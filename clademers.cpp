/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Dean M. Bobo
 *
 * CladeMers -- exact, lineage-specific k-mers from BLAST databases or FASTA.
 *
 * Read this file from top to bottom: options -> processes -> taxonomy -> k-mers ->
 * record parsing -> parallel subtraction -> outputs -> main.
 *
 * The biological operation is deliberately simple:
 *     result = k-mers in target references - k-mers in background references
 * A target is the requested taxid OR any descendant in the supplied taxonomy. Unknown
 * and mixed target/background annotations are background. Reverse complements
 * are equivalent. "Specific" means specific within THIS database snapshot;
 * a marker need not occur in every member of the target lineage.
 *
 * Build: c++ -O3 -std=c++17 -pthread clademers.cpp -lz -o clademers
 * BLAST sequence input uses blastdbcmd (NCBI BLAST+ 2.15 or newer).
 * All-plain FASTA/map/dump inputs need no external runtime program.
 * Indexed gzip FASTA uses zlib; compressed maps/dumps use the gzip CLI.
 * The optional BLAST taxonomy mode also uses the sqlite3 command-line reader.
 * Requires zlib headers/library; no Python or 128-bit extension is needed.
 * POSIX process handling targets Linux and macOS; use WSL on Windows.
 */

#include <algorithm>
#include <array>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <deque>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <memory>
#include <mutex>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// These are operating-system headers, not extra packages to install.
#include <fcntl.h>
#include <signal.h>
#include <spawn.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
extern char** environ;

namespace fs = std::filesystem;
using Taxid = std::uint32_t;
constexpr const char* VERSION = "0.13.0";
constexpr std::size_t CHUNK_BASES = 1024 * 1024;

volatile sig_atomic_t stop_requested = 0;
void request_stop(int) { stop_requested = 1; }
struct PauseRequested {}; // Deliberately not a runtime_error: parsing wrappers must preserve it.
void check_stop() { if (stop_requested) throw PauseRequested{}; }

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

std::string_view trim(std::string_view text) {
    const auto first = text.find_first_not_of(" \t\r\n");
    if (first == text.npos) return {};
    return text.substr(first, text.find_last_not_of(" \t\r\n") - first + 1);
}

std::uint64_t number(std::string_view text, const std::string& label) {
    if (text.empty()) fail("Missing integer: " + label);
    std::uint64_t value = 0;
    for (const char c : text) {
        if (c < '0' || c > '9') fail("Invalid integer for " + label + ": " + std::string(text));
        const unsigned digit = static_cast<unsigned>(c - '0');
        if (value > (UINT64_MAX - digit) / 10) fail("Integer overflow: " + label);
        value = value * 10 + digit;
    }
    return value;
}

Taxid taxid_number(std::string_view text) {
    const auto value = number(trim(text), "taxid");
    if (value == 0 || value > UINT32_MAX) fail("Taxids must be in 1..4294967295");
    return static_cast<Taxid>(value);
}

// ---------- 1. Command-line interface ----------

struct Options {
    fs::path nodes, names, blast_tax, fasta, taxid_map, output, results;
    std::string database, blastdbcmd = "blastdbcmd", sqlite3 = "sqlite3", gzip = "gzip";
    Taxid target = 0;
    unsigned k = 42, threads = 1, progress_seconds = 30, checkpoint_seconds = 300;
    bool resume = false;
    bool direct_resume = true;
    bool report_descendants = false;
    bool write_separate_fastas = false;
    bool gzip_separate_fastas = false;
    std::unordered_set<std::string> supplied;
};

void help() {
    std::cout << R"(CladeMers: exact lineage-specific k-mer extraction

Usage:
  clademers -wsf-gz -results /path/to/completed/results -out compressed_fastas
  clademers -wsf -rd -results /path/to/completed/results -out descendant_fastas
  clademers -rd -results /path/to/completed/results -out descendant_report
  clademers -resume /path/to/existing/result_directory -t 30
  clademers -taxnodes nodes.dmp -taxnames names.dmp -taxid 10226 \
           -k 42 -t 30 -blastdb /path/to/nt -out placozoa_k42
  clademers -blasttax /path/to/blast-taxonomy -taxid 10226 \
            -k 42 -t 30 -blastdb /path/to/nt -out placozoa_k42
  clademers -fasta references.fasta.gz -taxidmap sequence_taxids.tsv.gz \
            -taxnodes nodes.dmp.gz -taxnames names.dmp.gz -taxid 10226 -out placozoa_k42

Choose exactly one sequence input:
  -blastdb PREFIX  Complete nucleotide BLAST database prefix, not one volume
  OR
  -fasta FILE      Plain or gzip-compressed nucleotide FASTA, read twice
                  Gzip is detected from file contents; regular files only
  -taxidmap FILE   Optional plain/gzip two-column ID-to-taxid(s) map for FASTA
                  Without a map, use taxid=, taxids=, or kraken:taxid| headers

Choose exactly one taxonomy input:
  -taxnodes FILE    NCBI nodes.dmp or nodes.dmp.gz, matching the reference snapshot
  -taxnames FILE    NCBI names.dmp or names.dmp.gz (scientific-name labels)
  OR
  -blasttax DIR     Directory containing taxonomy4blast.sqlite3, taxdb.bti,
                   and taxdb.btd (requires the sqlite3 command-line reader)

Required in either mode:
  -taxid INTEGER   Target lineage; descendants are included automatically

Optional:
  -k INTEGER       Exact k-mer length, 1..64 [42]
  -t INTEGER       Background scanning worker threads [1]
  -rd, --report-descendants  Write descendant_summary.tsv (no value needed)
  -wsf, --write-separate-fastas  Write one FASTA per assigned LCA taxid, plus index
  -wsf-gz, --write-separate-fastas-gz  Same, gzip-compressed; implies -wsf
  -results DIRECTORY  With -rd, -wsf or -wsf-gz: use saved results without NT
                      Uses markers.fasta, target_taxonomy.tsv, manifest.tsv,
                      and SUCCESS; default output: DIRECTORY/descendant_report
  -out DIRECTORY   New output directory [clademers_<taxid>_k<k>]
  -blastdbcmd PATH Path to the BLAST+ executable [blastdbcmd in PATH]
  -sqlite3 PATH    SQLite CLI for -blasttax mode [sqlite3 in PATH]
  -gzip PATH       Reader for compressed maps/dumps or replay-mode FASTA [gzip]
  -progress-seconds INTEGER  Seconds between progress updates, 1..3600 [30]
  -checkpoint-seconds INTEGER  Save interval, 1..86400; 0 disables saving [300]
  -resume-mode MODE   direct [default] or replay; chosen when starting a new run
  -resume DIRECTORY  Resume a checkpoint; original inputs/options are restored
                     May override -t, reader paths, progress/checkpoint intervals
  -h, --help       Show this help
  --version        Show the CladeMers version

Outputs: markers.fasta, markers.tsv, target_taxids.txt, manifest.tsv, SUCCESS.
         All runs save LCA taxids in FASTA headers and specificity_taxid/name
         in TSV, plus target_taxonomy.tsv for later reporting with -rd -results.
         -wsf adds taxon_fastas/taxid_<id>.fasta and taxon_fastas.tsv; only
         direct assignments are included, and zero-marker taxa have no file.
         -wsf-gz writes .fasta.gz files with zlib; no external gzip is needed.
         BLAST input also writes database_info.txt and blastdbcmd.log.
         Subprocesses are recorded in commands.log; SQLite writes taxonomy.log.
         External gzip writes gzip.log; no decompressed file is created.

Only A/C/G/T windows are used; reverse complements are equivalent. Unknown
or mixed taxonomy is background. Results are exclusive within the supplied
reference database, not proof of universal biological specificity.
Descendant labels use the LCA of all observed taxids, not just the first source;
they do not establish presence in every genome within the assigned clade.
Candidate k-mers and an optional FASTA taxid map are held in RAM. Input reading
and target loading are serial; -t parallelizes background k-mer scanning.
Progress goes to stderr: percent, elapsed time, and estimated time remaining
for the named operation, not the whole run. Unknown totals show ETA unavailable.
Direct resume seeks BLAST OIDs/plain FASTA offsets or indexed gzip restart
points. Replay mode retains sequential readers and full prefix verification.
SIGUSR1/SIGTERM/SIGINT
request a checkpoint and exit 75; SIGKILL requires the previous periodic save.
)";
}

Options parse_options(int argc, char** argv) {
    Options options;
    std::unordered_set<std::string> seen;
    for (int i = 1; i < argc; ++i) {
        std::string flag = argv[i];
        if (flag == "-rd") flag = "--report-descendants";
        if (flag == "-wsf") flag = "--write-separate-fastas";
        if (flag == "-wsf-gz") flag = "--write-separate-fastas-gz";
        if (!seen.insert(flag).second) fail("Repeated option: " + flag);
        if (flag == "--report-descendants") { options.report_descendants = true; continue; }
        if (flag == "--write-separate-fastas") { options.write_separate_fastas = true; continue; }
        if (flag == "--write-separate-fastas-gz") {
            options.write_separate_fastas = true; options.gzip_separate_fastas = true; continue;
        }
        if (i + 1 == argc) fail("Missing value for " + flag);
        const std::string value = argv[++i];
        if (value.empty()) fail("Empty value for " + flag);
        if (flag == "-taxnodes") options.nodes = value;
        else if (flag == "-taxnames") options.names = value;
        else if (flag == "-blasttax") options.blast_tax = value;
        else if (flag == "-sqlite3") options.sqlite3 = value;
        else if (flag == "-gzip") options.gzip = value;
        else if (flag == "-taxid") options.target = taxid_number(value);
        else if (flag == "-blastdb") options.database = value;
        else if (flag == "-fasta") options.fasta = value;
        else if (flag == "-taxidmap") options.taxid_map = value;
        else if (flag == "-out") options.output = value;
        else if (flag == "-results") options.results = value;
        else if (flag == "-resume") { options.output = value; options.resume = true; }
        else if (flag == "-resume-mode") {
            if (value != "direct" && value != "replay") fail("-resume-mode must be direct or replay");
            options.direct_resume = value == "direct";
        }
        else if (flag == "-blastdbcmd") options.blastdbcmd = value;
        else if (flag == "-k") {
            const auto k = number(value, flag);
            if (k < 1 || k > 64) fail("-k must be between 1 and 64");
            options.k = static_cast<unsigned>(k);
        } else if (flag == "-t") {
            const auto threads = number(value, flag);
            if (threads < 1 || threads > 1024) fail("-t must be between 1 and 1024");
            options.threads = static_cast<unsigned>(threads);
        } else if (flag == "-progress-seconds") {
            const auto seconds = number(value, flag);
            if (seconds < 1 || seconds > 3600) fail("-progress-seconds must be between 1 and 3600");
            options.progress_seconds = static_cast<unsigned>(seconds);
        } else if (flag == "-checkpoint-seconds") {
            const auto seconds = number(value, flag);
            if (seconds > 86400) fail("-checkpoint-seconds must be between 0 and 86400");
            options.checkpoint_seconds = static_cast<unsigned>(seconds);
        } else fail("Unknown option: " + flag + " (see --help)");
    }
    options.supplied = seen;
    if (!options.results.empty()) {
        if (!options.report_descendants && !options.write_separate_fastas)
            fail("-results requires -rd, -wsf or -wsf-gz (or their long aliases)");
        for (const auto& flag : seen)
            if (flag != "-results" && flag != "--report-descendants" && flag != "--write-separate-fastas" &&
                flag != "--write-separate-fastas-gz" && flag != "-out" && flag != "-progress-seconds")
                fail("With -results, use only -rd, -wsf, -wsf-gz, -out and -progress-seconds; omit " + flag);
        options.results = fs::absolute(options.results);
        if (options.output.empty()) options.output = options.results / "descendant_report";
        options.output = fs::absolute(options.output);
        return options;
    }
    if (options.resume) {
        for (const auto& flag : seen)
            if (flag != "-resume" && flag != "-t" && flag != "-progress-seconds" &&
                flag != "-checkpoint-seconds" && flag != "-blastdbcmd" && flag != "-sqlite3" && flag != "-gzip" &&
                flag != "--report-descendants" && flag != "--write-separate-fastas" && flag != "--write-separate-fastas-gz")
                fail("With -resume, inputs, -taxid, -k and -out come from the checkpoint; omit " + flag);
        options.output = fs::absolute(options.output);
        return options;
    }
    if (!options.target) fail("Required: -taxid (see --help)");
    if (options.database.empty() == options.fasta.empty())
        fail("Choose exactly one sequence input: -blastdb PREFIX OR -fasta FILE");
    if (!options.taxid_map.empty() && options.fasta.empty())
        fail("-taxidmap is only used with -fasta");
    if (!options.fasta.empty() && seen.count("-blastdbcmd"))
        fail("-blastdbcmd is only used with -blastdb");
    if (options.fasta.empty() && !options.blast_tax.empty() && seen.count("-gzip"))
        fail("-gzip requires FASTA or .dmp taxonomy inputs");
    if (!options.blast_tax.empty()) {
        if (!options.nodes.empty() || !options.names.empty())
            fail("Choose either -blasttax DIR or -taxnodes FILE -taxnames FILE, not both");
    } else {
        if (options.nodes.empty() || options.names.empty())
            fail("Required taxonomy: -blasttax DIR OR both -taxnodes FILE and -taxnames FILE");
        if (seen.count("-sqlite3")) fail("-sqlite3 is only used with -blasttax");
    }
    if (options.database.find_first_of("\"\r\n") != std::string::npos)
        fail("BLAST database prefixes containing quotes or newlines are not supported");
    if (options.output.empty()) options.output = "clademers_" + std::to_string(options.target)
                                                + "_k" + std::to_string(options.k);
    // A resumed job may start in a different working directory.
    for (auto* path : {&options.nodes, &options.names, &options.blast_tax, &options.fasta, &options.taxid_map, &options.output})
        if (!path->empty()) *path = fs::absolute(*path);
    if (options.database.find('/') != std::string::npos) options.database = fs::absolute(options.database).string();
    for (auto* reader : {&options.blastdbcmd, &options.sqlite3, &options.gzip})
        if (reader->find('/') != std::string::npos) *reader = fs::absolute(*reader).string();
    return options;
}

// Progress is local to one operation. A total is optional: never infer a gzip
// stream's full size from its 32-bit footer, or treat four unequal stages as 25% each.
// Updates occur at work boundaries; no reporting thread or per-k-mer clock check.
std::string duration_text(double seconds) {
    if (!std::isfinite(seconds) || seconds < 0 || seconds >= static_cast<double>(UINT64_MAX))
        return "unavailable";
    const auto whole = static_cast<std::uint64_t>(seconds);
    std::ostringstream out;
    if (whole >= 86400) out << whole / 86400 << 'd' << ' ';
    out << std::setfill('0') << std::setw(2) << (whole / 3600) % 24 << ':'
        << std::setw(2) << (whole / 60) % 60 << ':' << std::setw(2) << whole % 60;
    return out.str();
}

class Progress {
    using Clock = std::chrono::steady_clock;
    std::string label_, unit_;
    std::optional<std::uint64_t> total_;
    std::uint64_t done_ = 0, initial_ = 0;
    unsigned interval_;
    Clock::time_point started_ = Clock::now(), last_ = started_;

    void report(Clock::time_point now, bool complete) const {
        const double elapsed = std::chrono::duration<double>(now - started_).count();
        std::ostringstream out;
        out << "      [progress] " << label_ << ": ";
        if (complete) out << "100.0%";
        else if (total_) out << std::fixed << std::setprecision(1)
                            << (total_.value() ? std::min(99.9, 100.0 * done_ / *total_) : 0.0) << '%';
        else out << "percent unavailable";
        out << " | elapsed " << duration_text(elapsed) << " | ETA ";
        if (complete) out << "00:00:00";
        else if (!total_) out << "unavailable (total unknown)";
        else if (done_ >= *total_) out << "finishing";
        else if (done_ <= initial_ || elapsed < 1.0) out << "estimating";
        else out << '~' << duration_text(std::ceil(elapsed * static_cast<double>(*total_ - done_) / (done_ - initial_)));
        out << " | " << done_;
        if (complete) out << '/' << done_;
        else if (total_) out << '/' << *total_;
        out << ' ' << unit_ << '\n';
        std::cerr << out.str();
    }
public:
    Progress(std::string label, std::string unit, std::optional<std::uint64_t> total, unsigned interval)
        : label_(std::move(label)), unit_(std::move(unit)), total_(total), interval_(interval) {
        report(started_, false);
    }
    void unknown_total() { total_.reset(); }
    void resume_at(std::uint64_t done) {
        done_ = initial_ = done;
        started_ = last_ = Clock::now();
        report(started_, false);
    }
    void update(std::uint64_t done) {
        done_ = done;
        // Changed inputs or approximate metadata must not yield >100% or negative ETA.
        if (total_ && done_ > *total_) total_.reset();
        const auto now = Clock::now();
        if (now - last_ >= std::chrono::seconds(interval_)) {
            report(now, false);
            last_ = now;
        }
    }
    void finish(std::uint64_t done) {
        done_ = done;
        report(Clock::now(), true); // Call only after successful reader/worker completion.
    }
};

// ---------- 2. Checked subprocesses, without a shell ----------

// We pass an argument vector directly to posix_spawnp. A database path may
// contain spaces or shell metacharacters; none are interpreted as commands.
// Reading EOF is NOT sufficient: finish() also checks the process exit code.
class Process {
    pid_t pid_ = -1;
    int fd_ = -1;
    std::array<char, 65536> buffer_{};
    std::size_t begin_ = 0, end_ = 0;
    std::string program_, log_name_;
public:
    Process(const std::vector<std::string>& arguments, const fs::path& output,
            const std::string& log_name = "blastdbcmd.log")
        : program_(arguments.at(0)), log_name_(log_name) {
        std::ofstream commands(output / "commands.log", std::ios::app);
        for (const auto& argument : arguments) commands << std::quoted(argument) << ' ';
        commands << '\n';
        commands.close();
        if (!commands) fail("Cannot write commands.log");
        int pipes[2];
        if (pipe(pipes) < 0) fail("Cannot create subprocess output pipe");
        int log = open((output / log_name_).c_str(), O_WRONLY | O_CREAT | O_APPEND, 0666);
        if (log < 0) { close(pipes[0]); close(pipes[1]); fail("Cannot open " + log_name_); }
        posix_spawn_file_actions_t actions;
        int error = posix_spawn_file_actions_init(&actions);
        if (error) {
            close(pipes[0]); close(pipes[1]); close(log);
            fail("Cannot initialize subprocess: " + std::string(std::strerror(error)));
        }
        auto action = [&](int result) { if (!error) error = result; };
        action(posix_spawn_file_actions_adddup2(&actions, pipes[1], STDOUT_FILENO));
        action(posix_spawn_file_actions_adddup2(&actions, log, STDERR_FILENO));
        action(posix_spawn_file_actions_addclose(&actions, pipes[0]));
        action(posix_spawn_file_actions_addclose(&actions, pipes[1]));
        action(posix_spawn_file_actions_addclose(&actions, log));
        std::vector<char*> argv;
        for (const auto& argument : arguments) argv.push_back(const_cast<char*>(argument.c_str()));
        argv.push_back(nullptr);
        if (!error) error = posix_spawnp(&pid_, argv[0], &actions, nullptr, argv.data(), environ);
        posix_spawn_file_actions_destroy(&actions);
        close(pipes[1]); close(log);
        if (error) {
            close(pipes[0]); pid_ = -1;
            fail("Cannot start " + arguments[0] + ": " + std::strerror(error)
                 + ". Check the executable path and PATH");
        }
        fd_ = pipes[0];
    }
    Process(const Process&) = delete;
    Process& operator=(const Process&) = delete;
    ~Process() {
        if (fd_ >= 0) close(fd_);
        if (pid_ > 0) {
            kill(pid_, SIGTERM); // Stop export if parsing/scanning failed.
            while (waitpid(pid_, nullptr, 0) < 0 && errno == EINTR) {}
        }
    }
    // BLAST/SQLite exports require complete newline-terminated rows. FASTA
    // permits a final unterminated line and needs its exact newline status.
    bool getline(std::string& line, bool* had_newline = nullptr) {
        check_stop();
        line.clear();
        if (had_newline) *had_newline = false;
        for (;;) {
            if (begin_ == end_) {
                ssize_t count;
                do {
                    count = read(fd_, buffer_.data(), buffer_.size());
                    if (count < 0 && errno == EINTR) check_stop();
                } while (count < 0 && errno == EINTR);
                check_stop();
                if (count < 0) fail("Error reading " + program_ + " output");
                if (count == 0) {
                    if (had_newline) return !line.empty();
                    if (!line.empty()) fail("Truncated " + program_ + " output: missing final newline");
                    return false;
                }
                begin_ = 0; end_ = static_cast<std::size_t>(count);
            }
            std::string_view available(buffer_.data() + begin_, end_ - begin_);
            const auto newline = available.find('\n');
            if (newline != available.npos) {
                line.append(available.data(), newline);
                begin_ += newline + 1;
                if (had_newline) *had_newline = true;
                return true;
            }
            line.append(available.data(), available.size());
            begin_ = end_;
        }
    }
    void finish() {
        check_stop();
        if (fd_ >= 0) { close(fd_); fd_ = -1; }
        int status = 0;
        pid_t result;
        do {
            result = waitpid(pid_, &status, 0);
            if (result < 0 && errno == EINTR) check_stop();
        } while (result < 0 && errno == EINTR);
        if (result < 0) fail("Cannot obtain " + program_ + " exit status");
        pid_ = -1;
        if (!WIFEXITED(status) || WEXITSTATUS(status) != 0)
            fail(program_ + " did not complete successfully; see " + log_name_ + " (no results accepted)");
    }
};

std::string capture(const std::vector<std::string>& arguments, const fs::path& output,
                    const std::string& log_name = "blastdbcmd.log") {
    Process process(arguments, output, log_name);
    std::string line, result;
    while (process.getline(line)) {
        result += line + '\n';
        if (result.size() > 16 * 1024 * 1024) fail("Unexpectedly large subprocess metadata output");
    }
    process.finish();
    if (result.empty()) fail(arguments.at(0) + " returned empty metadata");
    return result;
}

// ---------- 3. Taxonomy and small provenance helpers ----------

// FNV-1a is a small, non-cryptographic fingerprint for exact input bytes
// (decompressed bytes for gzip text inputs). It is for accidental-change detection,
// not authentication or identity proofs. Sequence identity NEVER depends on it.
struct Fingerprint {
    std::uint64_t value = 14695981039346656037ULL, bytes = 0;
    void add(std::string_view text) {
        for (unsigned char c : text) { value ^= c; value *= 1099511628211ULL; }
        bytes += text.size();
    }
    std::string hex() const {
        std::ostringstream out;
        out << std::hex << std::setw(16) << std::setfill('0') << value;
        return out.str();
    }
};

struct GzipPoint {
    std::uint64_t compressed = 0, uncompressed = 0, member_bytes = 0, crc = 0, bits = 0;
    std::string dictionary;
};

struct FastaCursor {
    std::uint64_t line = 0, records = 0, unknown_records = 0;
    std::string next_header; // Already fingerprinted lookahead; no completed sequence is retained.
    std::shared_ptr<const GzipPoint> gzip;
};

// Portable gzip access points follow zlib's zran approach: a DEFLATE block
// boundary, pending input bits, and the last 32 KiB of output. Save the running
// member CRC/length too, so a resumed raw inflater still validates gzip trailers.
// No temporary decompressed FASTA or architecture-specific zlib state is saved.
class IndexedGzip {
    std::ifstream file_;
    z_stream stream_{};
    std::array<unsigned char, 65536> input_{}, output_{};
    std::size_t in_begin_ = 0, in_end_ = 0, out_begin_ = 0, out_end_ = 0;
    std::uint64_t input_end_offset_ = 0, total_out_ = 0, consumed_ = 0, member_bytes_ = 0;
    uLong crc_ = 0;
    bool initialized_ = false, member_finished_ = false, eof_ = false;
    std::string history_;
    std::shared_ptr<const GzipPoint> ready_, pending_;
    std::uint64_t last_point_out_ = 0;

    bool input_available() {
        if (in_begin_ != in_end_) return true;
        check_stop();
        file_.read(reinterpret_cast<char*>(input_.data()), input_.size());
        if (file_.bad()) fail("gzip input read error");
        in_begin_ = 0; in_end_ = static_cast<std::size_t>(file_.gcount());
        input_end_offset_ += in_end_;
        return in_end_ != 0;
    }
    unsigned byte() {
        if (!input_available()) fail("Truncated gzip stream");
        return input_[in_begin_++];
    }
    std::uint32_t little32() {
        std::uint32_t value = 0;
        for (unsigned i = 0; i < 4; ++i) value |= static_cast<std::uint32_t>(byte()) << (8 * i);
        return value;
    }
    std::uint64_t input_position() const { return input_end_offset_ - (in_end_ - in_begin_); }
    void select_point() {
        if (pending_ && pending_->uncompressed <= consumed_) { ready_ = pending_; pending_.reset(); }
    }
    void add_point(unsigned bits) {
        auto point = std::make_shared<GzipPoint>();
        point->compressed = input_position(); point->uncompressed = total_out_;
        point->member_bytes = member_bytes_; point->crc = crc_; point->bits = bits;
        point->dictionary = history_;
        last_point_out_ = total_out_;
        if (total_out_ <= consumed_) ready_ = std::move(point);
        else pending_ = std::move(point);
    }
    bool start_member(bool first) {
        if (!input_available()) {
            if (first) fail("Empty gzip stream");
            return false;
        }
        if (!first && input_[in_begin_] == 0) {
            // gzip permits zero padding after its final member.
            while (input_available()) if (byte() != 0) fail("Nonzero garbage after gzip padding");
            return false;
        }
        uLong header_crc = crc32(0L, Z_NULL, 0);
        auto header_byte = [&] {
            const unsigned char value = static_cast<unsigned char>(byte());
            header_crc = crc32(header_crc, &value, 1);
            return static_cast<unsigned>(value);
        };
        if (header_byte() != 0x1f || header_byte() != 0x8b || header_byte() != 8)
            fail("Invalid gzip member header");
        const unsigned flags = header_byte();
        if (flags & 0xe0) fail("Reserved flags in gzip header");
        for (unsigned i = 0; i < 6; ++i) header_byte(); // MTIME, XFL, OS.
        if (flags & 4) {
            const unsigned low = header_byte(), high = header_byte();
            for (unsigned i = 0; i < low + (high << 8); ++i) header_byte();
        }
        if (flags & 8) while (header_byte()) {}
        if (flags & 16) while (header_byte()) {}
        if (flags & 2) {
            const unsigned low = byte(), high = byte();
            if ((header_crc & 0xffff) != low + (high << 8)) fail("gzip header CRC mismatch");
        }
        if (inflateReset(&stream_) != Z_OK) fail("Cannot reset gzip inflater");
        member_bytes_ = 0; crc_ = crc32(0L, Z_NULL, 0); history_.clear();
        member_finished_ = false;
        add_point(0);
        return true;
    }
    bool fill_output() {
        select_point();
        while (!eof_) {
            check_stop();
            if (member_finished_ && !start_member(false)) { eof_ = true; return false; }
            if (!input_available()) fail("Truncated gzip DEFLATE stream");
            stream_.next_in = input_.data() + in_begin_;
            stream_.avail_in = static_cast<uInt>(in_end_ - in_begin_);
            stream_.next_out = output_.data(); stream_.avail_out = output_.size();
            const auto before = stream_.avail_in;
            const int result = inflate(&stream_, Z_BLOCK);
            in_begin_ += before - stream_.avail_in;
            const auto produced = output_.size() - stream_.avail_out;
            if (result != Z_OK && result != Z_STREAM_END && result != Z_BUF_ERROR)
                fail("Invalid gzip DEFLATE data: " + std::string(stream_.msg ? stream_.msg : "inflate error"));
            // A pending end-of-block code can use only buffered bits: Z_BLOCK
            // may legitimately return at that boundary without consuming a byte.
            if (!produced && before == stream_.avail_in && result != Z_STREAM_END &&
                !(result == Z_OK && (stream_.data_type & 128)))
                fail("gzip inflater made no progress");
            crc_ = crc32(crc_, output_.data(), static_cast<uInt>(produced));
            member_bytes_ += produced; total_out_ += produced;
            if (produced >= 32768) history_.assign(reinterpret_cast<const char*>(output_.data() + produced - 32768), 32768);
            else if (produced) {
                if (history_.size() + produced > 32768) history_.erase(0, history_.size() + produced - 32768);
                history_.append(reinterpret_cast<const char*>(output_.data()), produced);
            }
            if (result == Z_STREAM_END) {
                const auto expected_crc = little32(), expected_size = little32();
                if (expected_crc != crc_ || expected_size != static_cast<std::uint32_t>(member_bytes_))
                    fail("gzip CRC or uncompressed-size mismatch");
                member_finished_ = true;
            } else if ((stream_.data_type & 128) && !(stream_.data_type & 64) && total_out_ - last_point_out_ >= CHUNK_BASES) {
                add_point(static_cast<unsigned>(stream_.data_type & 7));
            }
            out_begin_ = 0; out_end_ = produced;
            if (produced) return true;
        }
        return false;
    }
    void consume(std::size_t count) { out_begin_ += count; consumed_ += count; select_point(); }
public:
    explicit IndexedGzip(const fs::path& path) : file_(path, std::ios::binary) {
        if (!file_) fail("Cannot open indexed gzip input: " + path.string());
        if (inflateInit2(&stream_, -15) != Z_OK) fail("Cannot initialize gzip inflater");
        initialized_ = true;
        try { start_member(true); }
        catch (...) { inflateEnd(&stream_); initialized_ = false; throw; }
    }
    ~IndexedGzip() { if (initialized_) inflateEnd(&stream_); }
    IndexedGzip(const IndexedGzip&) = delete;
    std::shared_ptr<const GzipPoint> point() const { return ready_; }
    void seek(const std::shared_ptr<const GzipPoint>& point, std::uint64_t offset) {
        if (!point || point->bits > 7 || point->dictionary.size() > 32768 || point->crc > UINT32_MAX ||
            point->member_bytes < point->dictionary.size() || point->member_bytes > point->uncompressed ||
            point->uncompressed > offset || (point->bits && !point->compressed) ||
            point->compressed > static_cast<std::uint64_t>(std::numeric_limits<std::streamoff>::max()))
            fail("Invalid gzip checkpoint access point");
        file_.clear();
        file_.seekg(static_cast<std::streamoff>(point->compressed - (point->bits ? 1 : 0)));
        if (!file_) fail("Cannot seek gzip checkpoint position");
        in_begin_ = in_end_ = out_begin_ = out_end_ = 0;
        input_end_offset_ = point->compressed - (point->bits ? 1 : 0);
        if (inflateReset(&stream_) != Z_OK) fail("Cannot reset gzip inflater");
        if (point->bits && inflatePrime(&stream_, static_cast<int>(point->bits), byte() >> (8 - point->bits)) != Z_OK)
            fail("Cannot restore gzip checkpoint bits");
        if (!point->dictionary.empty() && inflateSetDictionary(&stream_,
            reinterpret_cast<const Bytef*>(point->dictionary.data()), point->dictionary.size()) != Z_OK)
            fail("Cannot restore gzip checkpoint dictionary");
        total_out_ = consumed_ = point->uncompressed; member_bytes_ = point->member_bytes;
        crc_ = static_cast<uLong>(point->crc); history_ = point->dictionary;
        ready_ = point; pending_.reset(); last_point_out_ = point->uncompressed;
        eof_ = member_finished_ = false;
        while (consumed_ < offset) {
            if (out_begin_ == out_end_ && !fill_output()) fail("gzip input ended before checkpoint offset");
            consume(static_cast<std::size_t>(std::min<std::uint64_t>(offset - consumed_, out_end_ - out_begin_)));
        }
        std::cerr << "[resume] Indexed gzip seek: compressed byte " << point->compressed
                  << "; replayed " << offset - point->uncompressed << " content bytes to the record boundary\n";
    }
    bool getline(std::string& line, bool& newline) {
        line.clear(); newline = false;
        for (;;) {
            if (out_begin_ == out_end_ && !fill_output()) return !line.empty();
            const std::string_view available(reinterpret_cast<const char*>(output_.data() + out_begin_), out_end_ - out_begin_);
            const auto end = available.find('\n');
            const auto count = end == available.npos ? available.size() : end;
            line.append(available.data(), count);
            consume(count + (end == available.npos ? 0 : 1));
            if (end != available.npos) { newline = true; return true; }
        }
    }
};

Fingerprint fingerprint_file(const fs::path& path, const Options& options) {
    std::ifstream file(path, std::ios::binary);
    if (!file) fail("Cannot open taxonomy file: " + path.string());
    Progress progress("Fingerprint " + path.filename().string(), "bytes", fs::file_size(path), options.progress_seconds);
    Fingerprint result;
    std::array<char, 65536> buffer{};
    while (file.read(buffer.data(), buffer.size()) || file.gcount()) {
        result.add(std::string_view(buffer.data(), static_cast<std::size_t>(file.gcount())));
        progress.update(result.bytes);
    }
    if (file.bad()) fail("Read error: " + path.string());
    progress.finish(result.bytes);
    return result;
}

// Regular files can be reopened for subsequent passes and input verification.
std::ifstream open_regular_input(const fs::path& path) {
    if (path == "-" || !fs::is_regular_file(path))
        fail("Expected a regular input file (not stdin or a pipe): " + path.string());
    std::ifstream file(path, std::ios::binary);
    if (!file) fail("Cannot open input file: " + path.string());
    return file;
}

bool is_gzip_input(std::ifstream& file, const fs::path& path) {
    unsigned char magic[2]{};
    file.read(reinterpret_cast<char*>(magic), 2);
    if (file.bad()) fail("Read error: " + path.string());
    const bool compressed = file.gcount() == 2 && magic[0] == 0x1f && magic[1] == 0x8b;
    file.clear();
    file.seekg(0);
    if (!file) fail("Cannot rewind input file: " + path.string());
    return compressed;
}

// Shared streaming reader for FASTA, maps, and dumps. Fingerprint content before any
// parsing/whitespace cleanup, and require a successful decompressor exit.
class TextInput {
    fs::path path_;
    std::ifstream file_;
    std::unique_ptr<Process> decompressor_;
    std::unique_ptr<IndexedGzip> indexed_;
    std::unique_ptr<Progress> progress_;
public:
    Fingerprint fingerprint;
    bool compressed;

    TextInput(const fs::path& path, const Options& options, const std::string& label,
              std::optional<std::uint64_t> known_bytes = std::nullopt, bool indexed = false)
        : path_(path), file_(open_regular_input(path)), compressed(is_gzip_input(file_, path)) {
        if (!compressed && path.extension() == ".gz")
            fail("File named .gz has no gzip header (invalid or truncated input): " + path.string());
        if (compressed) {
            file_.close();
            // Absolute filename after --: no shell, option injection, or output file.
            if (indexed) indexed_ = std::make_unique<IndexedGzip>(path);
            else decompressor_ = std::make_unique<Process>(
                    std::vector<std::string>{options.gzip, "-d", "-c", "--", fs::absolute(path).string()},
                    options.output, "gzip.log");
        }
        progress_ = std::make_unique<Progress>(label, "content bytes",
            known_bytes ? known_bytes : (compressed ? std::nullopt : std::optional<std::uint64_t>(fs::file_size(path))),
            options.progress_seconds);
    }

    bool getline(std::string& line) {
        check_stop();
        // The caller has finished parsing/processing the previously returned line.
        progress_->update(fingerprint.bytes);
        bool newline;
        if (indexed_) {
            if (!indexed_->getline(line, newline)) return false;
        } else if (decompressor_) {
            if (!decompressor_->getline(line, &newline)) return false;
        } else {
            const bool read_line = static_cast<bool>(std::getline(file_, line));
            check_stop(); // An interrupted file read must not look like the end of a FASTA record.
            if (!read_line) return false;
            newline = !file_.eof();
        }
        fingerprint.add(line);
        if (newline) fingerprint.add("\n");
        return true;
    }

    std::shared_ptr<const GzipPoint> access_point() const { return indexed_ ? indexed_->point() : nullptr; }
    void seek(const Fingerprint& prefix, const FastaCursor& cursor) {
        if (indexed_) indexed_->seek(cursor.gzip, prefix.bytes);
        else {
            if (compressed || prefix.bytes > static_cast<std::uint64_t>(std::numeric_limits<std::streamoff>::max()))
                fail("Input does not support direct checkpoint seeking");
            file_.clear(); file_.seekg(static_cast<std::streamoff>(prefix.bytes));
            if (!file_) fail("Cannot seek FASTA checkpoint offset");
            std::cerr << "[resume] Plain FASTA seek: byte " << prefix.bytes << '\n';
        }
        fingerprint = prefix;
        progress_->resume_at(prefix.bytes);
    }

    void finish() {
        try {
            if (decompressor_) decompressor_->finish();
            else if (file_.bad()) fail("Read error");
        } catch (const std::exception& e) {
            fail(path_.string() + ": " + e.what());
        }
        progress_->finish(fingerprint.bytes);
    }
};

struct TextFingerprint { Fingerprint fingerprint; bool compressed; };

TextFingerprint fingerprint_text(const fs::path& path, const Options& options, std::uint64_t known_bytes) {
    TextInput file(path, options, "Recheck " + path.filename().string(), known_bytes);
    std::string line;
    while (file.getline(line)) {} // Recheck bytes without loading another parsed data structure.
    file.finish();
    return {file.fingerprint, file.compressed};
}

// NCBI taxdb.bti has six big-endian uint32 header words, then sorted
// (taxid, byte offset) pairs into taxdb.btd. Its records contain four
// tab-separated name fields. Read bytes explicitly, independent of CPU endian.
// Format reference: NCBI seqdb_reader/seqdbtax.cpp (public-domain toolkit).
std::uint32_t read_be32(std::ifstream& file) {
    std::array<unsigned char, 4> bytes{};
    if (!file.read(reinterpret_cast<char*>(bytes.data()), bytes.size()))
        fail("Truncated taxdb.bti");
    std::uint32_t value = 0;
    for (const auto byte : bytes) value = (value << 8) | byte;
    return value;
}

std::unordered_map<Taxid, std::string> blast_taxon_names(
        const fs::path& directory, const std::unordered_set<Taxid>& requested) {
    const auto index_path = directory / "taxdb.bti", data_path = directory / "taxdb.btd";
    std::ifstream index(index_path, std::ios::binary), data(data_path, std::ios::binary);
    if (!index || !data) fail("Cannot open taxdb.bti and taxdb.btd in " + directory.string());
    if (read_be32(index) != 0x8739) fail("Unsupported or corrupt taxdb.bti: wrong magic number");
    const auto count = read_be32(index);
    for (unsigned i = 0; i < 4; ++i) read_be32(index); // Reserved header words.
    if (!count || fs::file_size(index_path) != 24ULL + 8ULL * count)
        fail("Corrupt taxdb.bti: entry count does not match file size");
    const auto data_size = fs::file_size(data_path);
    Taxid previous_id = 0;
    std::uint64_t previous_offset = 0;
    std::unordered_map<Taxid, std::string> names;
    auto read_name = [&](std::uint64_t end) {
        if (!requested.count(previous_id)) return;
        if (end <= previous_offset || end - previous_offset > 1024 * 1024)
            fail("Corrupt taxdb.btd: invalid target name-record length");
        data.seekg(static_cast<std::streamoff>(previous_offset));
        std::string record(static_cast<std::size_t>(end - previous_offset), '\0');
        if (!data.read(record.data(), static_cast<std::streamsize>(record.size())))
            fail("Truncated target name record in taxdb.btd");
        if (std::count(record.begin(), record.end(), '\t') != 3 ||
            record.find_first_of("\r\n") != record.npos || record.find('\0') != record.npos)
            fail("Malformed target name record in taxdb.btd");
        auto name = record.substr(0, record.find('\t'));
        if (trim(name).empty()) fail("Target has no scientific name in taxdb.btd");
        names.emplace(previous_id, std::move(name));
    };
    for (std::uint32_t i = 0; i < count; ++i) {
        const auto id = read_be32(index), offset = read_be32(index);
        if (!id || id <= previous_id || offset >= data_size || (i && offset <= previous_offset))
            fail("Corrupt taxdb.bti: unsorted taxids or invalid data offsets");
        if (i) read_name(offset);
        previous_id = id; previous_offset = offset;
    }
    // BLAST's name database legitimately omits some internal taxonomy nodes.
    // Membership comes from SQLite, so an absent display name is not an error.
    read_name(data_size);
    return names;
}

void check_sqlite_snapshot(const fs::path& path) {
    // Distributed taxonomy4blast files are standalone snapshots. Active SQLite
    // journals could make a fingerprint of the main file omit committed data.
    for (const auto* suffix : {"-wal", "-journal"}) {
        const fs::path sidecar = path.string() + suffix;
        if (fs::exists(sidecar) && fs::file_size(sidecar))
            fail("Use an inactive, standalone taxonomy4blast.sqlite3 snapshot (found "
                 + sidecar.string() + ")");
    }
}

// NCBI dump fields are separated by '|', with surrounding tabs/spaces.
std::vector<std::string_view> dump_fields(const std::string& line) {
    std::vector<std::string_view> fields;
    std::string_view rest(line);
    for (auto end = rest.find('|'); end != rest.npos; end = rest.find('|')) {
        fields.push_back(trim(rest.substr(0, end)));
        rest.remove_prefix(end + 1);
    }
    return fields;
}

template<class Function>
TextFingerprint read_dump(const fs::path& path, const Options& options, Function consume) {
    TextInput file(path, options, "Read " + path.filename().string());
    std::string line;
    std::uint64_t row = 0;
    while (file.getline(line)) {
        ++row;
        try { consume(dump_fields(line)); }
        catch (const std::exception& e) {
            fail(path.string() + ":" + std::to_string(row) + ": " + e.what());
        }
    }
    file.finish();
    return {file.fingerprint, file.compressed};
}

struct Taxonomy {
    struct Input {
        std::string key;
        fs::path path;
        Fingerprint fingerprint;
        bool text = false, compressed = false;
    };
    std::unordered_set<Taxid> targets;
    struct Clade { Taxid parent; std::size_t depth; std::string name; };
    // Retain the selected subtree to label every marker, even without -rd.
    std::unordered_map<Taxid, Clade> clades;
    std::string name, sqlite_version, name_source = "names.dmp";
    std::vector<Input> inputs;
    std::size_t node_count = 0;
    fs::path sqlite_path;

    Taxid lca(Taxid a, Taxid b) const {
        if (!a) return b; // No observations yet.
        if (!b) return a;
        while (a != b) {
            const auto& left = clades.at(a);
            const auto& right = clades.at(b);
            if (left.depth >= right.depth) a = left.parent;
            if (right.depth >= left.depth) b = right.parent;
        }
        return a;
    }

    void verify_unchanged(const Options& options) const {
        if (!sqlite_path.empty()) check_sqlite_snapshot(sqlite_path);
        for (const auto& input : inputs) {
            const auto current = input.text ? fingerprint_text(input.path, options, input.fingerprint.bytes)
                                            : TextFingerprint{fingerprint_file(input.path, options), false};
            if (current.compressed != input.compressed ||
                current.fingerprint.bytes != input.fingerprint.bytes ||
                current.fingerprint.value != input.fingerprint.value)
                fail("Taxonomy file changed during run: " + input.path.string());
        }
    }

    Taxonomy() = default; // Also used when loading a saved subtree for reporting.
    Taxonomy(const Options& options) {
        struct Node {
            Taxid id, parent;
            std::size_t parent_index = 0;
            unsigned char state = 0; // 0 = unseen, 1 = current path, 2 = resolved
            bool inside = false;
            std::size_t depth = 0;
        };
        std::vector<Node> nodes;
        if (options.blast_tax.empty()) {
            const auto content = read_dump(options.nodes, options, [&](const auto& fields) {
                if (fields.size() < 2) fail("Expected nodes.dmp fields separated by '|'");
                nodes.push_back({taxid_number(fields[0]), taxid_number(fields[1])});
            });
            inputs.push_back({"nodes", fs::absolute(options.nodes), content.fingerprint, true, content.compressed});
        } else {
            const auto directory = fs::absolute(options.blast_tax);
            if (!fs::is_directory(directory)) fail("-blasttax must name a directory: " + directory.string());
            sqlite_path = directory / "taxonomy4blast.sqlite3";
            check_sqlite_snapshot(sqlite_path);
            for (const auto& input : std::vector<std::pair<std::string, std::string>>{
                    {"taxonomy_sqlite", "taxonomy4blast.sqlite3"},
                    {"taxdb_bti", "taxdb.bti"}, {"taxdb_btd", "taxdb.btd"}}) {
                const auto path = directory / input.second;
                if (!fs::is_regular_file(path)) fail("Missing BLAST taxonomy file: " + path.string());
                inputs.push_back({input.first, path, fingerprint_file(path, options)});
            }
            sqlite_version = capture({options.sqlite3, "-init", "/dev/null", "-version"},
                                     options.output, "taxonomy.log");
            // Read the full tree and validate it using the same code as .dmp mode.
            // No recursive SQL/leaf-only filtering: internal nodes also matter.
            // -init prevents a user's .sqliterc from changing the export format.
            Progress progress("Read SQLite taxonomy", "nodes", std::nullopt, options.progress_seconds);
            Process reader({options.sqlite3, "-init", "/dev/null", "-readonly", "-batch", "-bail",
                "-noheader", "-separator", "\t", sqlite_path.string(),
                "PRAGMA quick_check; SELECT taxid, parent FROM TaxidInfo;"}, options.output, "taxonomy.log");
            std::string line;
            if (!reader.getline(line) || line != "ok")
                fail("SQLite taxonomy integrity check failed; see taxonomy.log");
            while (reader.getline(line)) {
                const auto tab = line.find('\t');
                if (tab == line.npos || line.find('\t', tab + 1) != line.npos)
                    fail("Malformed TaxidInfo export; expected taxid and parent");
                nodes.push_back({taxid_number(std::string_view(line).substr(0, tab)),
                                 taxid_number(std::string_view(line).substr(tab + 1))});
                progress.update(nodes.size());
            }
            reader.finish(); // Failed or truncated SQLite exports are never accepted.
            progress.finish(nodes.size());
        }
        std::sort(nodes.begin(), nodes.end(), [](const Node& a, const Node& b) { return a.id < b.id; });
        if (nodes.empty()) fail("Empty taxonomy tree");
        node_count = nodes.size();
        auto find_node = [&](Taxid id) {
            auto found = std::lower_bound(nodes.begin(), nodes.end(), id,
                [](const Node& node, Taxid key) { return node.id < key; });
            if (found == nodes.end() || found->id != id)
                fail("Taxid absent from taxonomy tree: " + std::to_string(id));
            return static_cast<std::size_t>(found - nodes.begin());
        };
        find_node(options.target); // A stale or mistyped target must not silently succeed.
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            auto& node = nodes[i];
            if (i && nodes[i - 1].id == node.id) fail("Duplicate taxid in taxonomy tree");
            node.parent_index = find_node(node.parent);
            if (node.id == node.parent) {
                if (node.id != 1) fail("Only NCBI root taxid 1 may be its own parent");
                node.state = 2;
                node.inside = (options.target == 1);
            }
        }
        // Follow parent links, then resolve the path from root toward children.
        // Memoization visits each link once after the parent-index lookup step.
        // Iteration avoids recursion limits for deep or malformed taxonomies.
        std::vector<std::size_t> path;
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            path.clear();
            std::size_t current = i;
            while (nodes[current].state == 0) {
                nodes[current].state = 1;
                path.push_back(current);
                current = nodes[current].parent_index;
            }
            if (nodes[current].state == 1) fail("Cycle in taxonomy tree");
            for (auto it = path.rbegin(); it != path.rend(); ++it) {
                auto& node = nodes[*it];
                node.inside = (node.id == options.target) || nodes[node.parent_index].inside;
                node.depth = nodes[node.parent_index].depth + 1;
                node.state = 2;
            }
        }
        for (const auto& node : nodes) if (node.inside) targets.insert(node.id);
        clades.reserve(targets.size());
        for (const auto& node : nodes) if (node.inside)
            clades.emplace(node.id, Clade{node.parent, node.depth, "taxid " + std::to_string(node.id)});
        std::unordered_map<Taxid, std::string> names;
        if (options.blast_tax.empty()) {
            const auto content = read_dump(options.names, options, [&](const auto& fields) {
                if (fields.size() < 4) fail("Expected four names.dmp fields");
                const Taxid id = taxid_number(fields[0]);
                if ((id == options.target || clades.count(id)) && fields[3] == "scientific name") {
                    if (!names.emplace(id, fields[1]).second) fail("Duplicate scientific name for target taxid");
                }
            });
            if (names[options.target].empty()) fail("Target has no scientific name in names.dmp");
            inputs.push_back({"names", fs::absolute(options.names), content.fingerprint, true, content.compressed});
        } else {
            names = blast_taxon_names(fs::absolute(options.blast_tax), targets);
            name_source = names.count(options.target) ? "taxdb" : "taxid_only";
        }
        name = names.count(options.target) ? names.at(options.target) : "taxid " + std::to_string(options.target);
        if (name_source == "taxid_only")
            std::cerr << "      No scientific name for this taxid in taxdb; using " << name << " as the label\n";
        for (const auto& entry : names) if (clades.count(entry.first) && !entry.second.empty()) {
            if (entry.second.find_first_of("\t\r\n") != std::string::npos || entry.second.find('\0') != std::string::npos)
                fail("Scientific name contains an invalid control character");
            clades.at(entry.first).name = entry.second;
        }
    }
};

// ---------- 4. Exact DNA representation (two ordinary 64-bit words) ----------

// A=00, C=01, G=10, T=11. At k=42 we need 84 bits; at k=64 we need 128.
// The full sequence is kept. Hash-table collisions are resolved by operator==.
struct Kmer {
    std::uint64_t high = 0, low = 0;
    bool operator==(const Kmer& other) const { return high == other.high && low == other.low; }
    bool operator<(const Kmer& other) const {
        return high < other.high || (high == other.high && low < other.low);
    }
};

std::uint64_t mix(std::uint64_t value) {
    value ^= value >> 30; value *= 0xbf58476d1ce4e5b9ULL;
    value ^= value >> 27; value *= 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

struct KmerHash {
    std::size_t operator()(const Kmer& key) const { return mix(key.low ^ mix(key.high)); }
};

std::uint64_t low_bits(unsigned bits) {
    return bits == 64 ? UINT64_MAX : (std::uint64_t(1) << bits) - 1;
}

int base_code(char base) {
    switch (base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

template<class Function>
void each_kmer(std::string_view sequence, unsigned k, Function consume) {
    const unsigned bits = 2 * k;
    const std::uint64_t mask_low = low_bits(std::min(bits, 64U));
    const std::uint64_t mask_high = bits > 64 ? low_bits(bits - 64) : 0;
    const unsigned complement_shift = bits - 2;
    Kmer forward, reverse;
    unsigned valid_bases = 0;
    for (std::size_t i = 0; i < sequence.size(); ++i) {
        const int base = base_code(sequence[i]);
        if (base < 0) {
            forward = reverse = {};
            valid_bases = 0; // No k-mer may bridge an ambiguous base.
            continue;
        }
        forward.high = ((forward.high << 2) | (forward.low >> 62)) & mask_high;
        forward.low = ((forward.low << 2) | static_cast<unsigned>(base)) & mask_low;
        reverse.low = (reverse.low >> 2) | (reverse.high << 62);
        reverse.high >>= 2;
        const auto complement = static_cast<std::uint64_t>(3 - base);
        if (complement_shift < 64) reverse.low |= complement << complement_shift;
        else reverse.high |= complement << (complement_shift - 64);
        if (valid_bases < k) ++valid_bases;
        if (valid_bases == k) {
            const bool reversed = reverse < forward;
            consume(reversed ? reverse : forward, i + 1 - k, reversed);
        }
    }
}

std::string decode(Kmer key, unsigned k) {
    std::string sequence(k, 'A');
    for (unsigned i = k; i > 0; --i) {
        sequence[i - 1] = "ACGT"[key.low & 3];
        key.low = (key.low >> 2) | (key.high << 62);
        key.high >>= 2;
    }
    return sequence;
}

// ---------- 5. Sequence input and classification using ALL reported taxonomy ----------

struct Record {
    std::string accession, taxids, leaf_taxids;
    std::string_view sequence; // View into the current line/record buffer; no copy.
    bool any_target = false, only_target = true;
    std::uint64_t input_record = 0; // 1-based FASTA record; zero for BLAST exports.
    std::uint64_t oid = UINT64_MAX;
    std::string_view export_line{};
    Taxid specificity_taxid = 0;
};

void classify_record(Record& record, const Taxonomy& taxonomy) {
    auto inspect_taxids = [&](std::string_view text) {
        bool found_any = false;
        while (!text.empty()) {
            const auto separator = text.find_first_of(";, \t");
            const auto token = trim(text.substr(0, separator));
            if (!token.empty()) {
                found_any = true;
                if (token == "N/A" || token == "0") record.only_target = false;
                else {
                    const Taxid id = taxid_number(token);
                    const bool inside = taxonomy.targets.count(id) != 0;
                    record.any_target |= inside;
                    record.only_target &= inside;
                    if (inside)
                        record.specificity_taxid = taxonomy.lca(record.specificity_taxid, id);
                }
            }
            if (separator == text.npos) break;
            text.remove_prefix(separator + 1);
        }
        if (!found_any) record.only_target = false;
    };
    inspect_taxids(record.taxids);
    // BLAST commonly reports N/A (or 0) when no extra leaf annotations exist.
    // That is acceptable IF the principal taxid is valid and inside the clade.
    const auto leaves = trim(record.leaf_taxids);
    if (!leaves.empty() && leaves != "N/A" && leaves != "0") inspect_taxids(leaves);
}

Record parse_record(std::string_view line, const Taxonomy& taxonomy) {
    std::array<std::string_view, 4> fields;
    std::string_view rest(line);
    for (unsigned i = 0; i < 3; ++i) {
        const auto tab = rest.find('\t');
        if (tab == rest.npos) fail("Malformed BLAST row: expected accession, taxid, leaf-taxids, sequence");
        fields[i] = rest.substr(0, tab);
        rest.remove_prefix(tab + 1);
    }
    if (rest.find('\t') != rest.npos) fail("Too many columns in BLAST row");
    fields[3] = rest;
    if (!fields[3].empty() && fields[3].back() == '\r') fields[3].remove_suffix(1);
    if (fields[0].empty() || fields[3].empty()) fail("Missing accession or sequence in BLAST row");
    Record record{std::string(fields[0]), std::string(fields[1]), std::string(fields[2]), fields[3]};
    classify_record(record, taxonomy);
    return record;
}

std::vector<std::string_view> words(std::string_view text) {
    std::vector<std::string_view> result;
    while (!(text = trim(text)).empty()) {
        const auto end = text.find_first_of(" \t\r\n");
        result.push_back(text.substr(0, end));
        if (end == text.npos) break;
        text.remove_prefix(end + 1);
    }
    return result;
}

using Taxids = std::vector<Taxid>;

void normalize_taxids(Taxids& ids) {
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
}

Taxids parse_taxid_list(std::string_view text) {
    Taxids ids;
    for (;;) {
        const auto separator = text.find_first_of(",;");
        const auto token = text.substr(0, separator);
        if (token.empty()) fail("Empty taxid in FASTA annotation or mapping list");
        ids.push_back(token == "N/A" || token == "0" ? 0 : taxid_number(token));
        if (separator == text.npos) break;
        text.remove_prefix(separator + 1);
    }
    normalize_taxids(ids);
    return ids;
}

std::string taxid_list_text(const Taxids& ids) {
    std::string text;
    for (const auto id : ids) {
        if (!text.empty()) text += ';';
        text += std::to_string(id);
    }
    return text;
}

struct FastaHeader { std::string id; Taxids taxids; };

FastaHeader parse_fasta_header(std::string_view header) {
    for (const unsigned char c : header)
        if ((c < 32 && c != '\t') || c == 127)
            fail("Control character in FASTA header; merged deflines must be converted to explicit taxids lists");
    const auto fields = words(header);
    if (fields.empty()) fail("FASTA header has no sequence ID");
    FastaHeader result{std::string(fields[0]), {}};
    auto accept = [&](Taxids ids) {
        if (!result.taxids.empty() && result.taxids != ids)
            fail("Conflicting taxid annotations in FASTA header: " + result.id);
        result.taxids = std::move(ids);
    };
    const std::string_view tag = "kraken:taxid|";
    for (auto pos = fields[0].find(tag); pos != fields[0].npos;
         pos = fields[0].find(tag, pos + tag.size())) {
        if (pos && fields[0][pos - 1] != '|') continue;
        const auto begin = pos + tag.size(), end = fields[0].find('|', begin);
        const auto value = fields[0].substr(begin, end == fields[0].npos ? end : end - begin);
        if (value.find_first_of(",;") != value.npos) fail("Kraken taxid annotations require one taxid");
        accept(parse_taxid_list(value));
    }
    for (std::size_t i = 1; i < fields.size(); ++i) {
        auto token = fields[i];
        const bool bracket = token.front() == '[';
        if (bracket) token.remove_prefix(1);
        const auto equals = token.find('=');
        if (equals == token.npos || (token.substr(0, equals) != "taxid" && token.substr(0, equals) != "taxids"))
            continue;
        if (bracket) {
            if (token.back() != ']') fail("Unclosed taxid annotation in FASTA header: " + result.id);
            token.remove_suffix(1);
        }
        accept(parse_taxid_list(token.substr(equals + 1)));
    }
    return result;
}

struct FastaTaxidMap {
    fs::path path;
    Fingerprint fingerprint;
    bool compressed;
    std::unordered_map<std::string, Taxids> entries;

    explicit FastaTaxidMap(const Options& options) : path(fs::absolute(options.taxid_map)) {
        TextInput file(options.taxid_map, options, "Read taxid map");
        compressed = file.compressed;
        std::string line;
        std::uint64_t row = 0;
        bool first_data = true;
        while (file.getline(line)) {
            ++row;
            const auto text = trim(line);
            if (text.empty() || text.front() == '#') continue;
            try {
                for (const unsigned char c : text)
                    if ((c < 32 && c != '\t') || c == 127) fail("Control character in taxid map");
                const auto fields = words(text);
                if (fields.size() != 2) fail("Expected two columns: sequence_id and taxid(s)");
                if (first_data && fields[0] == "sequence_id" &&
                    (fields[1] == "taxid" || fields[1] == "taxids")) { first_data = false; continue; }
                first_data = false;
                const auto ids = parse_taxid_list(fields[1]);
                auto& assigned = entries[std::string(fields[0])];
                assigned.insert(assigned.end(), ids.begin(), ids.end());
                normalize_taxids(assigned); // Repeated IDs preserve all assignments.
            } catch (const std::exception& e) {
                fail(path.string() + ":" + std::to_string(row) + ": " + e.what());
            }
        }
        file.finish();
        fingerprint = file.fingerprint;
        if (entries.empty()) fail("Empty FASTA taxid mapping file: " + path.string());
    }

    void verify_unchanged(const Options& options) const {
        const auto current = fingerprint_text(path, options, fingerprint.bytes);
        if (current.compressed != compressed || current.fingerprint.bytes != fingerprint.bytes ||
            current.fingerprint.value != fingerprint.value)
            fail("FASTA taxid mapping file changed during run: " + path.string());
    }
};

struct FastaScan {
    Fingerprint fingerprint;
    std::uint64_t records = 0, unknown_records = 0;
    bool compressed = false;
};

template<class Function>
FastaScan read_fasta(const Options& options, const Taxonomy& taxonomy,
                     const FastaTaxidMap* mapping, Function consume, const std::string& label,
                     std::optional<std::uint64_t> known_bytes = std::nullopt,
                     const FastaCursor* saved = nullptr, const Fingerprint* saved_prefix = nullptr) {
    const auto& path = options.fasta;
    TextInput file(path, options, label, known_bytes, options.direct_resume);
    FastaScan stats;
    stats.compressed = file.compressed;
    std::string line, sequence;
    FastaHeader header;
    std::uint64_t line_number = 0;
    if (saved) {
        if (!saved_prefix || saved->records == 0) fail("Invalid FASTA checkpoint cursor");
        file.seek(*saved_prefix, *saved);
        stats.records = saved->records; stats.unknown_records = saved->unknown_records; line_number = saved->line;
        if (!saved->next_header.empty()) header = parse_fasta_header(std::string_view(saved->next_header).substr(1));
    }
    auto emit = [&](const std::string& next_header) {
        if (header.id.empty()) return;
        if (sequence.empty()) fail("Empty FASTA sequence: " + header.id);
        Taxids ids = header.taxids;
        if (mapping) {
            const auto found = mapping->entries.find(header.id);
            if (found == mapping->entries.end()) ids = {0}; // Map mode never silently falls back to headers.
            else {
                if (!ids.empty() && ids != found->second)
                    fail("Header and mapping taxids disagree for FASTA ID: " + header.id);
                ids = found->second;
            }
        }
        if (ids.empty()) ids = {0};
        if (std::find(ids.begin(), ids.end(), 0) != ids.end()) ++stats.unknown_records;
        Record record{header.id, taxid_list_text(ids), "N/A", sequence};
        record.input_record = ++stats.records;
        classify_record(record, taxonomy);
        const FastaCursor cursor{line_number, stats.records, stats.unknown_records, next_header, file.access_point()};
        consume(record, file.fingerprint, cursor);
    };
    try {
        while (file.getline(line)) {
            ++line_number;
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (trim(line).empty()) continue;
            if (line.front() == '>') {
                emit(line);
                header = parse_fasta_header(std::string_view(line).substr(1));
                sequence.clear();
            } else {
                if (header.id.empty()) fail("Expected a FASTA header before sequence data");
                for (const unsigned char c : line) {
                    if (c == ' ' || c == '\t') continue;
                    const char upper = c >= 'a' && c <= 'z' ? static_cast<char>(c - 'a' + 'A') : static_cast<char>(c);
                    if (std::string_view("ACGTURYSWKMBDHVN-").find(upper) == std::string_view::npos)
                        fail("Invalid nucleotide character in FASTA sequence: " + header.id);
                    sequence.push_back(static_cast<char>(c)); // Ambiguities/gaps interrupt k-mer windows.
                }
            }
        }
        emit("");
        if (!stats.records) fail("Empty FASTA input");
        file.finish(); // CRC/truncation/exit errors invalidate all output.
        stats.fingerprint = file.fingerprint;
    } catch (const std::exception& e) {
        fail(path.string() + ":" + std::to_string(line_number) + ": " + e.what());
    }
    return stats;
}

struct Source { std::string accession, taxids, leaf_taxids; std::uint64_t input_record = 0; };
struct Origin { std::size_t source, position, index; bool reverse; Taxid specificity_taxid = 0; };
using Candidates = std::unordered_map<Kmer, Origin, KmerHash>;

// ---------- 6. Parallel, exact background subtraction ----------

// This immutable two-bit-check prefilter can say "definitely absent" or
// "check the hash table". It can NEVER remove a candidate by itself.
class Prefilter {
    std::vector<std::uint64_t> words_;
    std::size_t mask_;
public:
    explicit Prefilter(const Candidates& candidates) {
        std::size_t bits = 1024;
        if (candidates.size() > std::numeric_limits<std::size_t>::max() / 32)
            fail("Too many candidate k-mers for this platform");
        while (bits < candidates.size() * 16) bits *= 2;
        mask_ = bits - 1;
        words_.resize(bits / 64, 0);
        for (const auto& candidate : candidates) {
            const auto hash = KmerHash{}(candidate.first);
            for (const auto h : {hash, static_cast<std::size_t>(mix(hash))}) {
                const auto bit = h & mask_;
                words_[bit / 64] |= std::uint64_t(1) << (bit % 64);
            }
        }
    }
    bool maybe_contains(const Kmer& key) const {
        const auto hash = KmerHash{}(key);
        for (const auto h : {hash, static_cast<std::size_t>(mix(hash))}) {
            const auto bit = h & mask_;
            if (!(words_[bit / 64] & (std::uint64_t(1) << (bit % 64)))) return false;
        }
        return true;
    }
};

// The candidate table is read-only during scanning. Workers only set atomic
// exclusion flags, so simultaneous hits never race or change table structure.
// The bounded queue prevents workers from buffering the entire NT export.
class Scanner {
    const Candidates& candidates_;
    const Prefilter filter_;
    std::unique_ptr<std::atomic<bool>[]> excluded_;
    unsigned k_;
    std::size_t capacity_;
    std::deque<std::string> queue_;
    std::vector<std::thread> workers_;
    std::mutex mutex_;
    std::condition_variable ready_, space_, idle_;
    std::size_t active_ = 0;
    bool done_ = false, cancelled_ = false;
    std::exception_ptr error_;

    void worker() {
        try {
            for (;;) {
                std::string chunk;
                {
                    std::unique_lock<std::mutex> lock(mutex_);
                    ready_.wait(lock, [&] { return cancelled_ || done_ || !queue_.empty(); });
                    if (cancelled_ || (done_ && queue_.empty())) return;
                    chunk = std::move(queue_.front()); queue_.pop_front();
                    ++active_;
                    space_.notify_one();
                }
                each_kmer(chunk, k_, [&](const Kmer& key, std::size_t, bool) {
                    if (!filter_.maybe_contains(key)) return;
                    const auto found = candidates_.find(key); // All 2*k bits compared here.
                    if (found != candidates_.end())
                        excluded_[found->second.index].store(true, std::memory_order_relaxed);
                });
                {
                    std::lock_guard<std::mutex> lock(mutex_);
                    --active_;
                    if (!active_ && queue_.empty()) idle_.notify_all();
                }
            }
        } catch (...) {
            std::lock_guard<std::mutex> lock(mutex_);
            if (!error_) error_ = std::current_exception();
            cancelled_ = true;
            ready_.notify_all(); space_.notify_all(); idle_.notify_all();
        }
    }
public:
    Scanner(const Candidates& candidates, unsigned k, unsigned threads,
            const std::vector<unsigned char>& saved_flags = {})
        : candidates_(candidates), filter_(candidates),
          excluded_(new std::atomic<bool>[candidates.size()]), k_(k), capacity_(2 * threads) {
        if (!saved_flags.empty() && saved_flags.size() != (candidates.size() + 7) / 8)
            fail("Checkpoint exclusion flag count does not match candidates");
        for (std::size_t i = 0; i < candidates.size(); ++i)
            excluded_[i].store(!saved_flags.empty() && (saved_flags[i / 8] & (1U << (i % 8))));
        try {
            for (unsigned i = 0; i < threads; ++i) workers_.emplace_back([this] { worker(); });
        } catch (...) { cancel(); throw; }
    }
    Scanner(const Scanner&) = delete;
    Scanner& operator=(const Scanner&) = delete;
    ~Scanner() { cancel(); }
    void cancel() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            cancelled_ = true;
            ready_.notify_all(); space_.notify_all();
        }
        for (auto& thread : workers_) if (thread.joinable()) thread.join();
    }
    void add(std::string_view sequence) {
        // Each chunk includes k-1 bases of right overlap. Therefore every
        // window crossing a chunk boundary is scanned, once, in the left chunk.
        for (std::size_t start = 0; start + k_ <= sequence.size(); start += CHUNK_BASES) {
            std::string chunk(sequence.substr(start, CHUNK_BASES + k_ - 1));
            std::unique_lock<std::mutex> lock(mutex_);
            space_.wait(lock, [&] { return cancelled_ || queue_.size() < capacity_; });
            if (error_) std::rethrow_exception(error_);
            if (cancelled_) fail("Background workers cancelled");
            queue_.push_back(std::move(chunk));
            ready_.notify_one();
        }
    }
    void finish() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            done_ = true; ready_.notify_all();
        }
        for (auto& thread : workers_) if (thread.joinable()) thread.join();
        if (error_) std::rethrow_exception(error_);
    }
    void drain() {
        std::unique_lock<std::mutex> lock(mutex_);
        idle_.wait(lock, [&] { return error_ || (queue_.empty() && !active_); });
        if (error_) std::rethrow_exception(error_);
    }
    bool excluded(std::size_t index) const { return excluded_[index].load(std::memory_order_relaxed); }
};

// ---------- 7. Auditable output and failure handling ----------

void write_text(const fs::path& path, const std::string& text) {
    std::ofstream out(path, std::ios::binary);
    if (!out) fail("Cannot open output: " + path.string());
    out << text;
    out.close();
    if (!out) fail("Output write failed (disk full?): " + path.string());
}

std::string utc_now() {
    const auto now = std::time(nullptr);
    std::tm utc{};
    gmtime_r(&now, &utc);
    std::ostringstream result;
    result << std::put_time(&utc, "%Y-%m-%dT%H:%M:%SZ");
    return result.str();
}

// A two-column manifest stays readable and needs no JSON dependency. Escape
// literal tabs/newlines/backslashes in values so every entry remains one row.
std::string tsv_escape(std::string_view value) {
    std::string result;
    for (char c : value) {
        if (c == '\\') result += "\\\\";
        else if (c == '\t') result += "\\t";
        else if (c == '\n') result += "\\n";
        else if (c == '\r') result += "\\r";
        else result += c;
    }
    return result;
}

class Manifest {
    fs::path directory_;
    std::vector<std::pair<std::string, std::string>> values_;
public:
    explicit Manifest(fs::path directory) : directory_(std::move(directory)) {}
    const auto& entries() const { return values_; }
    std::string get(const std::string& key) const {
        for (const auto& entry : values_) if (entry.first == key) return entry.second;
        return {};
    }
    void set(const std::string& key, const std::string& value) {
        for (auto& entry : values_) if (entry.first == key) { entry.second = value; return; }
        values_.emplace_back(key, value);
    }
    void count(const std::string& key, std::uint64_t value) { set(key, std::to_string(value)); }
    void save() const {
        std::string text = "key\tvalue\n";
        for (const auto& entry : values_) text += entry.first + '\t' + tsv_escape(entry.second) + '\n';
        write_text(directory_ / "manifest.partial.tsv", text);
        fs::rename(directory_ / "manifest.partial.tsv", directory_ / "manifest.tsv");
    }
};

using TaxonCounts = std::unordered_map<Taxid, std::uint64_t>;

std::vector<Taxid> sorted_taxids(const Taxonomy& taxonomy) {
    std::vector<Taxid> ids(taxonomy.targets.begin(), taxonomy.targets.end());
    std::sort(ids.begin(), ids.end());
    return ids;
}

Fingerprint write_taxonomy_snapshot(const fs::path& directory, const Taxonomy& taxonomy) {
    std::ofstream out(directory / "target_taxonomy.partial.tsv", std::ios::binary);
    if (!out) fail("Cannot open target taxonomy snapshot");
    Fingerprint fingerprint;
    auto line = [&](const std::string& text) { out << text; fingerprint.add(text); };
    line("taxid\tparent_taxid\tscientific_name\n");
    for (Taxid id : sorted_taxids(taxonomy)) {
        const auto& clade = taxonomy.clades.at(id);
        line(std::to_string(id) + '\t' + std::to_string(clade.parent) + '\t' + clade.name + '\n');
        check_stop();
    }
    out.close();
    if (!out) fail("Taxonomy snapshot write failed (disk full?)");
    fs::rename(directory / "target_taxonomy.partial.tsv", directory / "target_taxonomy.tsv");
    return fingerprint;
}

void write_descendant_summary(const fs::path& directory, Taxid target,
                              const Taxonomy& taxonomy, const TaxonCounts& assigned) {
    auto inclusive = assigned;
    const auto ids = sorted_taxids(taxonomy);
    auto deepest_first = ids;
    std::sort(deepest_first.begin(), deepest_first.end(), [&](Taxid a, Taxid b) {
        return taxonomy.clades.at(a).depth > taxonomy.clades.at(b).depth;
    });
    for (Taxid id : deepest_first) {
        if (id != target) inclusive[taxonomy.clades.at(id).parent] += inclusive[id];
        check_stop();
    }
    std::ofstream summary(directory / "descendant_summary.partial.tsv");
    if (!summary) fail("Cannot open descendant summary");
    summary << "taxid\tparent_taxid\tscientific_name\tassigned_kmers\tassigned_kmers_incl_descendants\n";
    for (Taxid id : ids) {
        const auto& clade = taxonomy.clades.at(id);
        const auto count = assigned.find(id);
        summary << id << '\t' << clade.parent << '\t' << clade.name << '\t'
                << (count == assigned.end() ? 0 : count->second) << '\t' << inclusive[id] << '\n';
        check_stop();
    }
    summary.close();
    if (!summary) fail("Descendant summary write failed (disk full?)");
    check_stop();
    fs::rename(directory / "descendant_summary.partial.tsv", directory / "descendant_summary.tsv");
}

// Keep the number of open files bounded even for targets with thousands of
// descendant taxa. Reopened files append within this fresh staging directory.
// Publication happens only after every record and source fingerprint is valid.
class TaxonFastaWriter {
    static constexpr std::size_t MAX_OPEN = 32;
    struct Handle {
        std::ofstream file;
        gzFile gzip = nullptr;
        std::list<Taxid>::iterator recent;
        // On failure/interrupt, close any outstanding gzip members. Their
        // staging files are never published; a retry regenerates all outputs.
        ~Handle() { if (gzip) gzclose(gzip); }
    };
    fs::path directory_, staging_;
    const Taxonomy& taxonomy_;
    bool compressed_;
    TaxonCounts counts_;
    std::unordered_map<Taxid, Handle> open_;
    std::list<Taxid> recent_;

    std::string filename(Taxid id) const {
        return "taxid_" + std::to_string(id) + (compressed_ ? ".fasta.gz" : ".fasta");
    }

    static void clear_generated_directory(const fs::path& path) {
        const auto status = fs::symlink_status(path);
        if (!fs::exists(status)) return;
        if (!fs::is_directory(status)) fail("Expected a real generated FASTA directory: " + path.string());
        // Output-phase resume may replace an incomplete earlier export. Refuse
        // unexpected contents instead of recursively deleting unrelated files.
        std::vector<fs::path> files;
        for (const auto& entry : fs::directory_iterator(path)) {
            const auto name = entry.path().filename().string();
            // Recognize both formats when a resumed run enables compression.
            const auto plain = name.size() > 3 && name.compare(name.size() - 3, 3, ".gz") == 0
                ? name.substr(0, name.size() - 3) : name;
            if (!fs::is_regular_file(entry.symlink_status()) || plain.size() <= 12 ||
                plain.compare(0, 6, "taxid_") || plain.compare(plain.size() - 6, 6, ".fasta"))
                fail("Unexpected file in generated FASTA directory: " + entry.path().string());
            const auto id = taxid_number(std::string_view(plain).substr(6, plain.size() - 12));
            if (plain != "taxid_" + std::to_string(id) + ".fasta") fail("Unexpected FASTA filename: " + name);
            files.push_back(entry.path());
        }
        for (const auto& file : files) fs::remove(file);
        fs::remove(path);
    }

    void close_oldest() {
        const Taxid id = recent_.back();
        auto& handle = open_.at(id);
        if (handle.gzip) {
            auto file = handle.gzip; handle.gzip = nullptr;
            if (gzclose(file) != Z_OK) fail("Compressed FASTA close failed (disk full?): " + filename(id));
        } else {
            handle.file.close();
            if (!handle.file) fail("Separate FASTA write failed (disk full?): " + filename(id));
        }
        open_.erase(id); recent_.pop_back();
    }
public:
    TaxonFastaWriter(const fs::path& directory, const Taxonomy& taxonomy, bool compressed)
        : directory_(directory), staging_(directory / "taxon_fastas.partial"), taxonomy_(taxonomy), compressed_(compressed) {
        clear_generated_directory(staging_);
        if (!fs::create_directory(staging_)) fail("Cannot create separate FASTA staging directory");
    }

    void write(Taxid id, const std::string& record) {
        if (!taxonomy_.clades.count(id)) fail("Separate FASTA taxid is outside target subtree");
        auto found = open_.find(id);
        if (found == open_.end()) {
            if (open_.size() == MAX_OPEN) close_oldest();
            found = open_.try_emplace(id).first;
            auto& handle = found->second;
            if (compressed_) {
                // Eviction closes a complete member; reopening appends another
                // standard gzip member, preserving records without re-reading.
                handle.gzip = gzopen((staging_ / filename(id)).c_str(), "ab6");
                if (!handle.gzip) fail("Cannot open compressed FASTA: " + filename(id));
            } else {
                handle.file.open(staging_ / filename(id), std::ios::binary | std::ios::app);
                if (!handle.file) fail("Cannot open separate FASTA: " + filename(id));
            }
            recent_.push_front(id); handle.recent = recent_.begin();
        } else recent_.splice(recent_.begin(), recent_, found->second.recent);
        auto& handle = found->second;
        if (compressed_) {
            if (record.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
                fail("FASTA record exceeds gzip write limit");
            if (gzwrite(handle.gzip, record.data(), static_cast<unsigned>(record.size())) != static_cast<int>(record.size()))
                fail("Compressed FASTA write failed (disk full?): " + filename(id));
        } else {
            handle.file << record;
            if (!handle.file) fail("Separate FASTA write failed (disk full?): " + filename(id));
        }
        ++counts_[id];
    }

    std::size_t file_count() const { return counts_.size(); }

    void finish() {
        while (!open_.empty()) close_oldest();
        std::ofstream index(directory_ / "taxon_fastas.partial.tsv");
        if (!index) fail("Cannot open separate FASTA index");
        index << "taxid\tscientific_name\tassigned_kmers\tfasta_file\n";
        for (Taxid id : sorted_taxids(taxonomy_)) {
            const auto found = counts_.find(id);
            if (found == counts_.end()) continue;
            index << id << '\t' << taxonomy_.clades.at(id).name << '\t' << found->second
                  << "\ttaxon_fastas/" << filename(id) << '\n';
            check_stop();
        }
        index.close();
        if (!index) fail("Separate FASTA index write failed (disk full?)");
        check_stop();
        clear_generated_directory(directory_ / "taxon_fastas");
        fs::rename(staging_, directory_ / "taxon_fastas");
        fs::rename(directory_ / "taxon_fastas.partial.tsv", directory_ / "taxon_fastas.tsv");
    }
};

// Reporting reads only the portable result files. It never resolves the saved
// reference paths, starts sequence readers, or loads candidate checkpoints.
class SavedResultReader {
    fs::path path_;
    std::ifstream file_;
public:
    Fingerprint fingerprint;
    explicit SavedResultReader(const fs::path& path) : path_(path), file_(path, std::ios::binary) {
        if (!file_) fail("Cannot read saved result: " + path.string());
    }
    bool getline(std::string& line) {
        check_stop();
        if (!std::getline(file_, line)) {
            if (!file_.eof() || file_.bad()) fail("Error reading saved result: " + path_.string());
            return false;
        }
        if (file_.eof()) fail("Truncated saved result (missing final newline): " + path_.string());
        if (line.find('\0') != line.npos) fail("Invalid NUL in saved result: " + path_.string());
        fingerprint.add(line); fingerprint.add("\n");
        return true;
    }
};

std::vector<std::string> tab_fields(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t start = 0;
    for (;;) {
        const auto end = line.find('\t', start);
        fields.push_back(line.substr(start, end == line.npos ? end : end - start));
        if (end == line.npos) return fields;
        start = end + 1;
    }
}

void check_saved_fingerprint(const Fingerprint& actual, const Manifest& manifest, const std::string& key) {
    if (manifest.get(key + "_fnv1a64") != actual.hex() ||
        number(manifest.get(key + "_bytes"), key + " bytes") != actual.bytes)
        fail("Saved result fingerprint mismatch: " + key + " (file changed or damaged)");
}

Taxonomy read_saved_taxonomy(const fs::path& path, Taxid target, const Manifest& manifest) {
    SavedResultReader reader(path);
    std::string line;
    if (!reader.getline(line) || line != "taxid\tparent_taxid\tscientific_name")
        fail("Unsupported target_taxonomy.tsv header");
    Taxonomy taxonomy;
    Taxid previous = 0;
    std::unordered_map<Taxid, std::vector<Taxid>> children;
    while (reader.getline(line)) {
        const auto fields = tab_fields(line);
        if (fields.size() != 3 || trim(fields[2]).empty() || fields[2].find('\r') != fields[2].npos)
            fail("Malformed saved taxonomy row");
        const auto id = taxid_number(fields[0]), parent = taxid_number(fields[1]);
        if (id <= previous) fail("Saved taxonomy taxids must be unique and sorted");
        previous = id;
        taxonomy.targets.insert(id);
        taxonomy.clades.emplace(id, Taxonomy::Clade{parent, 0, fields[2]});
        if (id != target) children[parent].push_back(id);
    }
    check_saved_fingerprint(reader.fingerprint, manifest, "target_taxonomy");
    if (!taxonomy.targets.count(target) || taxonomy.targets.size() != number(manifest.get("target_taxids"), "target taxids"))
        fail("Saved taxonomy is missing the target or has the wrong node count");
    const auto root_parent = taxonomy.clades.at(target).parent;
    if ((target == 1 && root_parent != 1) || (target != 1 && taxonomy.targets.count(root_parent)))
        fail("Invalid saved target parent");
    // All nodes must be reachable from the selected root. Disconnected nodes,
    // absent parents and cycles cannot masquerade as part of the target clade.
    std::vector<Taxid> queue{target};
    for (std::size_t i = 0; i < queue.size(); ++i) {
        const auto found = children.find(queue[i]);
        if (found != children.end()) for (Taxid child : found->second) {
            taxonomy.clades.at(child).depth = taxonomy.clades.at(queue[i]).depth + 1;
            queue.push_back(child);
        }
        check_stop();
    }
    if (queue.size() != taxonomy.targets.size()) fail("Saved taxonomy has missing parents or disconnected cycles");
    taxonomy.name = taxonomy.clades.at(target).name;
    return taxonomy;
}

void report_saved_results(const Options& options) {
    if (!fs::is_regular_file(options.results / "SUCCESS"))
        fail("-results requires a completed extraction with SUCCESS");
    Manifest source(options.results);
    SavedResultReader metadata(options.results / "manifest.tsv");
    std::string line;
    if (!metadata.getline(line) || line != "key\tvalue") fail("Malformed saved manifest");
    std::unordered_set<std::string> keys;
    while (metadata.getline(line)) {
        const auto fields = tab_fields(line);
        if (fields.size() != 2 || fields[0].empty() || !keys.insert(fields[0]).second)
            fail("Malformed or duplicate saved manifest entry");
        source.set(fields[0], fields[1]); // Fields used here are numeric or fixed tokens, never escaped paths.
    }
    if (source.get("status") != "complete") fail("Saved extraction is not complete");
    if (source.get("results_format") != "1" || source.get("lca_saved") != "yes" || source.get("mode") != "extraction")
        fail("Saved results lack the portable LCA catalog; use results from CladeMers 0.11.0 or later. Older runs need a new extraction");
    const auto target = taxid_number(source.get("target_taxid"));
    const auto k = number(source.get("k"), "saved k"), expected = number(source.get("retained_kmers"), "retained markers");
    if (!k || k > 64) fail("Invalid saved k-mer length");
    std::cerr << "[report] Reading saved taxonomy and marker LCA labels; no reference scan\n";
    const auto taxonomy = read_saved_taxonomy(options.results / "target_taxonomy.tsv", target, source);
    std::unique_ptr<TaxonFastaWriter> separate;
    if (options.write_separate_fastas) {
        if (!fs::create_directory(options.output)) fail("Output directory already exists: " + options.output.string());
        separate = std::make_unique<TaxonFastaWriter>(options.output, taxonomy, options.gzip_separate_fastas);
    }
    SavedResultReader markers(options.results / "markers.fasta");
    TaxonCounts assigned;
    std::uint64_t count = 0;
    std::string previous, sequence;
    const auto prefix = ">km_" + std::to_string(target) + "_k" + std::to_string(k) + "_";
    Progress progress("Read saved marker labels", "markers", expected, options.progress_seconds);
    while (markers.getline(line)) {
        if (count >= expected) fail("Saved marker count exceeds manifest");
        const auto header = prefix + std::to_string(count + 1) + " taxid=";
        if (line.compare(0, header.size(), header) != 0) fail("Malformed saved marker header or missing LCA taxid");
        const Taxid id = taxid_number(std::string_view(line).substr(header.size()));
        if (!taxonomy.targets.count(id)) fail("Saved marker LCA is outside the selected subtree");
        if (!markers.getline(sequence) || sequence.size() != k || sequence.find_first_not_of("ACGT") != sequence.npos ||
            (!previous.empty() && sequence <= previous)) fail("Saved markers must be distinct, sorted ACGT k-mers of the saved length");
        each_kmer(sequence, static_cast<unsigned>(k), [&](const Kmer&, std::size_t, bool reverse) {
            if (reverse) fail("Saved marker is not in canonical orientation");
        });
        previous = sequence;
        if (separate) separate->write(id, line + '\n' + sequence + '\n');
        ++assigned[id]; ++count;
        if (count % 65536 == 0) progress.update(count);
    }
    if (count != expected) fail("Saved marker count does not match manifest");
    check_saved_fingerprint(markers.fingerprint, source, "markers_fasta");
    progress.finish(count);
    if (!options.write_separate_fastas && !fs::create_directory(options.output))
        fail("Output directory already exists: " + options.output.string());
    Manifest report(options.output);
    report.set("clademers_version", VERSION);
    report.set("mode", "descendant_report");
    report.set("report_descendants", options.report_descendants ? "yes" : "no");
    report.set("write_separate_fastas", options.write_separate_fastas ? "yes" : "no");
    report.set("separate_fasta_compression", options.gzip_separate_fastas ? "gzip" : "none");
    report.set("zlib_version", zlibVersion());
    report.set("source_results", fs::absolute(options.results).string());
    report.set("source_clademers_version", source.get("clademers_version"));
    report.set("status", "running");
    report.set("started_utc", utc_now());
    report.set("reference_scan", "none; saved marker LCA labels and taxonomy snapshot only");
    report.count("target_taxid", target); report.set("target_name", taxonomy.name);
    report.count("k", k); report.count("retained_kmers", count);
    report.count("target_taxids", taxonomy.targets.size());
    for (const auto* key : {"markers_fasta_fnv1a64", "markers_fasta_bytes", "target_taxonomy_fnv1a64", "target_taxonomy_bytes"})
        report.set(key, source.get(key));
    report.save();
    if (options.report_descendants) write_descendant_summary(options.output, target, taxonomy, assigned);
    if (separate) {
        separate->finish();
        report.count("separate_fasta_files", separate->file_count());
        report.set("separate_fasta_assignment", "one file per exact LCA; ancestors do not duplicate markers");
    }
    report.set("status", "complete"); report.set("finished_utc", utc_now()); report.save();
    write_text(options.output / "SUCCESS.partial", "Descendant report completed from saved LCA assignments; no reference scan.\n");
    fs::rename(options.output / "SUCCESS.partial", options.output / "SUCCESS");
    std::cerr << "Reported " << count << " saved markers across " << taxonomy.targets.size()
              << " taxa. Results: " << options.output << '\n';
}

// ---------- Durable, portable checkpoints ----------

// One writer per result directory; the kernel releases the advisory lock on exit,
// including SIGKILL. O_CLOEXEC keeps readers from inheriting the descriptor.
class RunLock {
    int fd_ = -1;
public:
    explicit RunLock(const fs::path& directory) {
        fd_ = open((directory / "run.lock").c_str(), O_RDWR | O_CREAT | O_CLOEXEC, 0666);
        if (fd_ < 0) fail("Cannot open run.lock");
        struct flock lock{};
        lock.l_type = F_WRLCK; lock.l_whence = SEEK_SET;
        if (fcntl(fd_, F_SETLK, &lock) < 0) {
            close(fd_); fd_ = -1;
            fail("Result directory is locked by another run (or filesystem locking is unavailable)");
        }
    }
    ~RunLock() { if (fd_ >= 0) close(fd_); }
    RunLock(const RunLock&) = delete;
};

void sync_path(const fs::path& path) {
    const int fd = open(path.c_str(), O_RDONLY | O_CLOEXEC);
    if (fd < 0) fail("Cannot open checkpoint for sync: " + path.string());
    int result;
    do { result = fsync(fd); } while (result < 0 && errno == EINTR);
    const int error = errno;
    close(fd);
    if (result < 0) fail("Cannot sync checkpoint: " + path.string() + ": " + std::strerror(error));
}

// Explicit little-endian uint64 values and length-prefixed strings: no native
// structs, pointers, size_t layout, or CPU-specific serialization. The trailing
// FNV fingerprint detects accidental damage; it is not authentication.
class CheckpointWriter {
    fs::path path_, temporary_;
    std::ofstream file_;
    Fingerprint fingerprint_;
public:
    explicit CheckpointWriter(const fs::path& path)
        : path_(path), temporary_(path.string() + ".partial"), file_(temporary_, std::ios::binary | std::ios::trunc) {
        if (!file_) fail("Cannot create checkpoint: " + temporary_.string());
    }
    void bytes(const char* data, std::size_t count) {
        file_.write(data, static_cast<std::streamsize>(count));
        if (!file_) fail("Checkpoint write failed (disk full?): " + temporary_.string());
        fingerprint_.add(std::string_view(data, count));
    }
    void integer(std::uint64_t value) {
        char data[8];
        for (unsigned i = 0; i < 8; ++i) data[i] = static_cast<char>(value >> (8 * i));
        bytes(data, 8);
    }
    void string(const std::string& value) { integer(value.size()); bytes(value.data(), value.size()); }
    void fingerprint(const Fingerprint& value) { integer(value.value); integer(value.bytes); }
    std::string finish() {
        const auto saved = fingerprint_;
        integer(saved.value); integer(saved.bytes);
        file_.close();
        if (!file_) fail("Checkpoint close failed (disk full?): " + temporary_.string());
        sync_path(temporary_);
        fs::rename(temporary_, path_); // The old checkpoint remains intact until this atomic replacement.
        sync_path(path_.parent_path());
        return saved.hex() + ":" + std::to_string(saved.bytes);
    }
};

class CheckpointReader {
    std::ifstream file_;
    std::uint64_t remaining_;
    Fingerprint fingerprint_;
public:
    explicit CheckpointReader(const fs::path& path) : file_(open_regular_input(path)), remaining_(fs::file_size(path)) {
        if (remaining_ < 16) fail("Truncated checkpoint: " + path.string());
    }
    void bytes(char* data, std::size_t count) {
        if (count > remaining_ || !file_.read(data, static_cast<std::streamsize>(count)))
            fail("Truncated checkpoint payload");
        remaining_ -= count;
        fingerprint_.add(std::string_view(data, count));
    }
    std::uint64_t integer() {
        unsigned char data[8]; bytes(reinterpret_cast<char*>(data), 8);
        std::uint64_t value = 0;
        for (unsigned i = 0; i < 8; ++i) value |= static_cast<std::uint64_t>(data[i]) << (8 * i);
        return value;
    }
    std::size_t count(std::uint64_t minimum_bytes) {
        const auto value = integer();
        if (value > remaining_ / minimum_bytes || value > std::numeric_limits<std::size_t>::max())
            fail("Invalid checkpoint length/count");
        return static_cast<std::size_t>(value);
    }
    std::string string() {
        const auto size = count(1);
        std::string value(size, '\0'); bytes(value.data(), size); return value;
    }
    Fingerprint fingerprint() {
        Fingerprint value; value.value = integer(); value.bytes = integer(); return value;
    }
    std::string finish() {
        const auto saved = fingerprint_;
        const auto expected = fingerprint();
        if (remaining_ || expected.value != saved.value || expected.bytes != saved.bytes || file_.peek() != EOF)
            fail("Checkpoint checksum/length mismatch; checkpoint is corrupt");
        return saved.hex() + ":" + std::to_string(saved.bytes);
    }
};

void write_options(CheckpointWriter& out, const Options& options) {
    for (const auto& path : {options.nodes, options.names, options.blast_tax, options.fasta, options.taxid_map})
        out.string(path.string());
    for (const auto& value : {options.database, options.blastdbcmd, options.sqlite3, options.gzip}) out.string(value);
    out.integer(options.target); out.integer(options.k); out.integer(options.threads);
    out.integer(options.progress_seconds); out.integer(options.checkpoint_seconds);
    out.integer(options.direct_resume);
    out.integer(options.report_descendants);
    out.integer(options.write_separate_fastas);
    out.integer(options.gzip_separate_fastas);
}

void read_options(CheckpointReader& in, Options& options) {
    for (auto* path : {&options.nodes, &options.names, &options.blast_tax, &options.fasta, &options.taxid_map}) *path = in.string();
    for (auto* value : {&options.database, &options.blastdbcmd, &options.sqlite3, &options.gzip}) *value = in.string();
    const auto target = in.integer(), k = in.integer(), threads = in.integer();
    const auto progress = in.integer(), checkpoint = in.integer();
    if (!target || target > UINT32_MAX || !k || k > 64 || !threads || threads > 1024 ||
        !progress || progress > 3600 || checkpoint > 86400 || options.database.empty() == options.fasta.empty() ||
        (options.blast_tax.empty() && (options.nodes.empty() || options.names.empty()))) fail("Invalid checkpoint options");
    options.target = static_cast<Taxid>(target); options.k = static_cast<unsigned>(k);
    options.threads = static_cast<unsigned>(threads); options.progress_seconds = static_cast<unsigned>(progress);
    options.checkpoint_seconds = static_cast<unsigned>(checkpoint);
    const auto direct = in.integer();
    if (direct > 1) fail("Invalid checkpoint resume mode");
    options.direct_resume = direct != 0;
    const auto descendants = in.integer();
    if (descendants > 1) fail("Invalid checkpoint descendant-report option");
    options.report_descendants = descendants != 0;
    const auto separate = in.integer();
    if (separate > 1) fail("Invalid checkpoint separate-FASTA option");
    options.write_separate_fastas = separate != 0;
    const auto compressed_fastas = in.integer();
    if (compressed_fastas > 1 || (compressed_fastas && !options.write_separate_fastas))
        fail("Invalid checkpoint separate-FASTA compression option");
    options.gzip_separate_fastas = compressed_fastas != 0;
}

void write_candidates(CheckpointWriter& out, const Candidates& candidates, const std::vector<Source>& sources) {
    out.integer(sources.size());
    for (const auto& source : sources) {
        out.string(source.accession); out.string(source.taxids); out.string(source.leaf_taxids); out.integer(source.input_record);
    }
    out.integer(candidates.size());
    for (const auto& item : candidates) {
        out.integer(item.first.high); out.integer(item.first.low);
        out.integer(item.second.source); out.integer(item.second.position);
        out.integer(item.second.index); out.integer(item.second.reverse);
        out.integer(item.second.specificity_taxid);
    }
}

void read_candidates(CheckpointReader& in, Candidates& candidates, std::vector<Source>& sources, unsigned k) {
    const auto source_count = in.count(32);
    sources.reserve(source_count);
    for (std::size_t i = 0; i < source_count; ++i) {
        Source source;
        source.accession = in.string(); source.taxids = in.string(); source.leaf_taxids = in.string(); source.input_record = in.integer();
        if (source.accession.empty()) fail("Invalid checkpoint source");
        sources.push_back(std::move(source));
    }
    const auto count = in.count(56);
    candidates.reserve(count);
    std::vector<bool> indices(count, false);
    for (std::size_t i = 0; i < count; ++i) {
        Kmer key; key.high = in.integer(); key.low = in.integer();
        const auto source = in.integer(), position = in.integer(), index = in.integer(), reverse = in.integer();
        const auto specificity = in.integer();
        if (source >= sources.size() || position > std::numeric_limits<std::size_t>::max() ||
            index >= count || indices[index] || reverse > 1 || specificity > UINT32_MAX ||
            (key.low & ~low_bits(std::min(2 * k, 64U))) ||
            (key.high & ~(k > 32 ? low_bits(2 * k - 64) : 0))) fail("Invalid checkpoint candidate");
        indices[index] = true;
        if (!candidates.emplace(key, Origin{static_cast<std::size_t>(source), static_cast<std::size_t>(position),
            static_cast<std::size_t>(index), reverse != 0, static_cast<Taxid>(specificity)}).second)
            fail("Duplicate checkpoint candidate");
    }
}

struct RunState {
    // 0: loading targets; 1: scanning background; 2: inputs verified, writing output.
    std::uint64_t phase = 0, rows = 0, target_rows = 0, target_bases = 0;
    std::uint64_t background_rows = 0, background_bases = 0, background_windows = 0;
    Fingerprint prefix;
    FastaScan first_scan;
    std::string identity, base_hash;
    std::vector<unsigned char> flags;
    std::uint64_t oid = UINT64_MAX, oid_rows = 0, blast_bases = 0, oid_length = 0;
    Fingerprint oid_prefix;
    FastaCursor fasta_cursor;
};

void write_cursor(CheckpointWriter& out, const RunState& state) {
    out.integer(state.oid); out.integer(state.oid_rows); out.integer(state.blast_bases); out.integer(state.oid_length);
    out.fingerprint(state.oid_prefix);
    const auto& cursor = state.fasta_cursor;
    out.integer(cursor.line); out.integer(cursor.records); out.integer(cursor.unknown_records); out.string(cursor.next_header);
    out.integer(cursor.gzip ? 1 : 0);
    if (cursor.gzip) {
        const auto& point = *cursor.gzip;
        out.integer(point.compressed); out.integer(point.uncompressed); out.integer(point.member_bytes);
        out.integer(point.crc); out.integer(point.bits); out.string(point.dictionary);
    }
}

void read_cursor(CheckpointReader& in, RunState& state) {
    state.oid = in.integer(); state.oid_rows = in.integer(); state.blast_bases = in.integer(); state.oid_length = in.integer();
    state.oid_prefix = in.fingerprint();
    auto& cursor = state.fasta_cursor;
    cursor.line = in.integer(); cursor.records = in.integer(); cursor.unknown_records = in.integer(); cursor.next_header = in.string();
    const auto has_point = in.integer();
    if (has_point > 1 || state.oid_rows > state.rows || cursor.unknown_records > cursor.records ||
        (state.oid_rows && state.oid >= INT32_MAX) || (!cursor.next_header.empty() && cursor.next_header.front() != '>'))
        fail("Invalid checkpoint cursor");
    if (has_point) {
        auto point = std::make_shared<GzipPoint>();
        point->compressed = in.integer(); point->uncompressed = in.integer(); point->member_bytes = in.integer();
        point->crc = in.integer(); point->bits = in.integer(); point->dictionary = in.string();
        if (point->bits > 7 || point->crc > UINT32_MAX || point->dictionary.size() > 32768 ||
            point->member_bytes < point->dictionary.size() || point->member_bytes > point->uncompressed ||
            point->uncompressed > state.prefix.bytes) fail("Invalid indexed gzip checkpoint");
        cursor.gzip = std::move(point);
    }
}

std::string phase_name(std::uint64_t phase) {
    return phase == 0 ? "target loading" : (phase == 1 ? "background scan" : "output");
}

void write_checkpoint(const Options& options, RunState& state, const Candidates& candidates,
                      const std::vector<Source>& sources, Scanner* scanner, const Manifest& manifest) {
    const auto directory = options.output / "checkpoint";
    fs::create_directory(directory);
    std::cerr << "      [checkpoint] Saving " << phase_name(state.phase) << " after " << state.rows << " records\n";
    if (scanner) scanner->drain(); // Exclusion flags and input position must describe the same completed work.
    if (state.phase && state.base_hash.empty()) {
        CheckpointWriter base(directory / "candidates.bin");
        base.string("CladeMers candidates 2"); base.string(VERSION);
        write_candidates(base, candidates, sources);
        state.base_hash = base.finish(); // Immutable thereafter; do not rewrite large candidate tables per scan checkpoint.
    }
    CheckpointWriter out(directory / "state.bin");
    out.string("CladeMers checkpoint 5"); out.string(VERSION);
    write_options(out, options);
    out.integer(state.phase); out.integer(state.rows); out.integer(state.target_rows); out.integer(state.target_bases);
    out.integer(state.background_rows); out.integer(state.background_bases); out.integer(state.background_windows);
    out.fingerprint(state.prefix); out.fingerprint(state.first_scan.fingerprint);
    out.integer(state.first_scan.records); out.integer(state.first_scan.unknown_records); out.integer(state.first_scan.compressed);
    out.string(state.identity); out.string(state.base_hash);
    write_cursor(out, state);
    out.integer(manifest.entries().size());
    for (const auto& entry : manifest.entries()) { out.string(entry.first); out.string(entry.second); }
    if (!state.phase) write_candidates(out, candidates, sources);
    else {
        out.integer((candidates.size() + 7) / 8);
        for (std::size_t i = 0; i < candidates.size(); i += 8) {
            unsigned char flags = 0;
            for (std::size_t bit = 0; bit < 8 && i + bit < candidates.size(); ++bit)
                if (scanner ? scanner->excluded(i + bit) :
                    (!state.flags.empty() && (state.flags[(i + bit) / 8] & (1U << ((i + bit) % 8)))))
                    flags |= static_cast<unsigned char>(1U << bit);
            out.bytes(reinterpret_cast<const char*>(&flags), 1);
        }
    }
    out.finish();
    std::cerr << "      [checkpoint] Saved " << phase_name(state.phase) << " after " << state.rows << " records\n";
}

void read_checkpoint(Options& options, RunState& state, Candidates& candidates,
                     std::vector<Source>& sources, Manifest& manifest) {
    const auto requested = options;
    const auto directory = options.output / "checkpoint";
    CheckpointReader in(directory / "state.bin");
    if (in.string() != "CladeMers checkpoint 5" || in.string() != VERSION)
        fail("Unsupported checkpoint format/version; resume with the same CladeMers release");
    read_options(in, options);
    state.phase = in.integer(); state.rows = in.integer(); state.target_rows = in.integer(); state.target_bases = in.integer();
    state.background_rows = in.integer(); state.background_bases = in.integer(); state.background_windows = in.integer();
    state.prefix = in.fingerprint(); state.first_scan.fingerprint = in.fingerprint();
    state.first_scan.records = in.integer(); state.first_scan.unknown_records = in.integer();
    const auto compressed = in.integer(); state.first_scan.compressed = compressed != 0;
    if (state.phase > 2 || compressed > 1 || state.background_windows > state.background_rows ||
        (state.phase && state.background_rows > state.rows) || (!state.phase && state.target_rows > state.rows))
        fail("Invalid checkpoint state");
    state.identity = in.string(); state.base_hash = in.string();
    read_cursor(in, state);
    const auto entries = in.count(16);
    for (std::size_t i = 0; i < entries; ++i) {
        const auto key = in.string(); const auto value = in.string(); manifest.set(key, value);
    }
    if (!state.phase) read_candidates(in, candidates, sources, options.k);
    else {
        const auto count = in.count(1);
        state.flags.resize(count); in.bytes(reinterpret_cast<char*>(state.flags.data()), count);
    }
    in.finish();
    if (state.phase) {
        CheckpointReader base(directory / "candidates.bin");
        if (base.string() != "CladeMers candidates 2" || base.string() != VERSION) fail("Unsupported candidate checkpoint");
        read_candidates(base, candidates, sources, options.k);
        if (base.finish() != state.base_hash || state.flags.size() != (candidates.size() + 7) / 8 || candidates.empty())
            fail("Candidate checkpoint does not match saved scan state");
    }
    if (requested.supplied.count("-t")) options.threads = requested.threads;
    if (requested.supplied.count("-progress-seconds")) options.progress_seconds = requested.progress_seconds;
    if (requested.supplied.count("-checkpoint-seconds")) options.checkpoint_seconds = requested.checkpoint_seconds;
    if (requested.supplied.count("-blastdbcmd")) options.blastdbcmd = requested.blastdbcmd;
    if (requested.supplied.count("-sqlite3")) options.sqlite3 = requested.sqlite3;
    if (requested.supplied.count("-gzip")) options.gzip = requested.gzip;
    if (requested.supplied.count("--report-descendants")) options.report_descendants = true;
    if (requested.supplied.count("--write-separate-fastas")) options.write_separate_fastas = true;
    if (requested.supplied.count("--write-separate-fastas-gz")) {
        options.write_separate_fastas = true; options.gzip_separate_fastas = true;
    }
    for (auto* reader : {&options.blastdbcmd, &options.sqlite3, &options.gzip})
        if (reader->find('/') != std::string::npos) *reader = fs::absolute(*reader).string();
}

bool blast_prefix_exists(const fs::path& prefix) {
    return fs::is_regular_file(prefix.string() + ".nal") || fs::is_regular_file(prefix.string() + ".nin");
}

fs::path resolve_blast_prefix(const std::string& name, const fs::path& directory, const Options& options) {
    fs::path path(name);
    if (path.extension() == ".nal") path.replace_extension();
    auto found = [&](const fs::path& candidate) -> fs::path {
        if (!blast_prefix_exists(candidate)) return {};
        return fs::weakly_canonical(fs::absolute(candidate));
    };
    if (const auto candidate = found(path.is_absolute() ? path : directory / path); !candidate.empty()) return candidate;
    if (!path.is_absolute()) {
        const auto search = capture({options.blastdbcmd, "-show_blastdb_search_path"}, options.output);
        std::string_view rest = trim(search);
        while (!rest.empty()) {
            const auto colon = rest.find(':');
            const fs::path base(std::string(rest.substr(0, colon)));
            if (const auto candidate = found((base.empty() ? fs::current_path() : base) / path); !candidate.empty()) return candidate;
            if (colon == rest.npos) break;
            rest.remove_prefix(colon + 1);
        }
    }
    fail("Cannot resolve local BLAST prefix for direct resume: " + name +
         ". Supply its absolute prefix, or start with -resume-mode replay");
}

std::vector<std::string> alias_words(std::string_view text) {
    std::vector<std::string> values;
    while (!(text = trim(text)).empty()) {
        if (text.front() == '"') {
            text.remove_prefix(1);
            const auto end = text.find('"');
            if (end == text.npos) fail("Unclosed quoted path in BLAST alias");
            values.emplace_back(text.substr(0, end)); text.remove_prefix(end + 1);
            if (!text.empty() && text.front() != ' ' && text.front() != '\t') fail("Malformed BLAST alias path list");
        } else {
            const auto end = text.find_first_of(" \t\r\n");
            values.emplace_back(text.substr(0, end));
            if (end == text.npos) break;
            text.remove_prefix(end);
        }
    }
    return values;
}

// Freeze the complete local alias/volume/filter-file set, not just -info text.
// Stat checks avoid rereading sequence payloads. They detect ordinary edits,
// truncation and replacement; immutable references remain a required assumption.
std::vector<fs::path> blast_files(const Options& options) {
    std::vector<fs::path> pending{options.database}, files;
    std::unordered_set<std::string> seen;
    while (!pending.empty()) {
        const auto prefix = pending.back(); pending.pop_back();
        if (!seen.insert(prefix.string()).second) continue;
        if (seen.size() > 100000) fail("BLAST alias graph is unexpectedly large");
        const auto directory = prefix.parent_path();
        const auto stem = prefix.filename().string() + ".";
        for (const auto& entry : fs::directory_iterator(directory)) {
            const auto name = entry.path().filename().string();
            if (name.compare(0, stem.size(), stem) == 0) {
                const auto suffix = name.substr(stem.size());
                if (!suffix.empty() && suffix.front() == 'n' && suffix.find('.') == suffix.npos && entry.is_regular_file())
                    files.push_back(entry.path());
            }
        }
        const fs::path alias(prefix.string() + ".nal");
        if (!fs::is_regular_file(alias)) {
            if (!fs::is_regular_file(prefix.string() + ".nin")) fail("Missing BLAST volume index: " + prefix.string());
            continue;
        }
        if (fs::file_size(alias) > 16 * 1024 * 1024) fail("BLAST alias file is unexpectedly large");
        std::ifstream input(alias);
        std::string line;
        bool have_databases = false;
        while (std::getline(input, line)) {
            const auto text = trim(line);
            if (text.empty() || text.front() == '#') continue;
            const auto split = text.find_first_of(" \t");
            const auto key = text.substr(0, split);
            const auto value = split == text.npos ? std::string_view{} : trim(text.substr(split));
            if (key == "DBLIST") {
                for (const auto& name : alias_words(value)) {
                    pending.push_back(resolve_blast_prefix(name, directory, options)); have_databases = true;
                }
            } else if (key == "GILIST" || key == "TILIST" || key == "SEQIDLIST" || key == "OIDLIST" || key == "TAXIDLIST") {
                const auto paths = alias_words(value);
                if (paths.size() != 1) fail("Expected one filter filename in BLAST alias: " + alias.string());
                const fs::path path(paths.front());
                files.push_back(path.is_absolute() ? path : directory / path);
            } else if (key != "TITLE" && key != "NSEQ" && key != "LENGTH" && key != "MAXOID" &&
                       key != "MAXLEN" && key != "STATS_TOTLEN" && key != "STATS_NSEQ" &&
                       key != "FIRST_OID" && key != "LAST_OID" && key != "MEMB_BIT" && key != "OID_MASK_TYPE")
                fail("Unsupported BLAST alias field for direct resume: " + std::string(key) + "; use -resume-mode replay for this database");
        }
        if (input.bad() || !have_databases) fail("Invalid BLAST alias DBLIST: " + alias.string());
    }
    std::sort(files.begin(), files.end()); files.erase(std::unique(files.begin(), files.end()), files.end());
    if (files.empty()) fail("No local BLAST files found for direct resume");
    return files;
}

std::string file_snapshot(const std::vector<fs::path>& paths) {
    std::ostringstream out;
    for (const auto& path : paths) {
        struct stat info{};
        if (stat(path.c_str(), &info) < 0 || !S_ISREG(info.st_mode)) fail("Cannot stat reference file: " + path.string());
#if defined(__APPLE__)
        const auto modified = info.st_mtimespec, changed = info.st_ctimespec;
#else
        const auto modified = info.st_mtim, changed = info.st_ctim;
#endif
        out << std::quoted(fs::canonical(path).string()) << ' ' << info.st_size << ' '
            << modified.tv_sec << ':' << modified.tv_nsec << ' ' << changed.tv_sec << ':' << changed.tv_nsec << '\n';
    }
    return out.str();
}

std::string sequence_snapshot(const Options& options) {
    if (!options.direct_resume) return {};
    return file_snapshot(options.fasta.empty() ? blast_files(options) : std::vector<fs::path>{options.fasta});
}

std::string input_identity(const Options& options, const Taxonomy& taxonomy, const FastaTaxidMap* mapping,
                            const std::string& database_info, const std::string& blast_version) {
    std::ostringstream out;
    out << std::quoted(database_info) << ' ' << std::quoted(blast_version) << '\n';
    for (const auto& input : taxonomy.inputs)
        out << input.key << ' ' << input.fingerprint.hex() << ' ' << input.fingerprint.bytes << ' ' << input.compressed << '\n';
    if (mapping) out << mapping->fingerprint.hex() << ' ' << mapping->fingerprint.bytes << ' ' << mapping->compressed << '\n';
    if (!options.fasta.empty()) {
        auto file = open_regular_input(options.fasta);
        out << "fasta gzip " << is_gzip_input(file, options.fasta) << '\n';
    }
    // BLAST prefix resolution can depend on the submission environment. Require
    // the same environment when a bare database name is used, as well as -info.
    if (!options.database.empty() && options.database.find('/') == std::string::npos) {
        const char* blastdb = std::getenv("BLASTDB");
        out << std::quoted(fs::current_path().string()) << ' ' << std::quoted(blastdb ? blastdb : "") << '\n';
    }
    out << sequence_snapshot(options);
    return out.str();
}

std::vector<std::string> blast_arguments(const Options& options, const std::string& prefix = {}) {
    // BLAST itself parses -db as a whitespace-separated list, independently of
    // shell parsing. Inner quotes preserve ONE prefix, including any spaces.
    return {options.blastdbcmd, "-db", "\"" + (prefix.empty() ? options.database : prefix) + "\"", "-dbtype", "nucl"};
}

std::string blast_resume_alias(const Options& options, const RunState& state) {
    if (!options.direct_resume || !state.rows) return {};
    if (!state.oid_rows || state.oid >= INT32_MAX) fail("Invalid BLAST OID checkpoint");
    const auto prefix = fs::absolute(options.output / "resume_range");
    // NCBI aliases use one-based FIRST_OID and inclusive LAST_OID, while %o
    // exports zero-based OIDs. Keep the original alias as DBLIST to retain ALL
    // of its filters and volume ordering. Revalidate only the saved OID group.
    write_text(prefix.string() + ".nal", "TITLE CladeMers resume range\nDBLIST \"" + options.database +
        "\"\nFIRST_OID " + std::to_string(state.oid + 1) + "\nLAST_OID 2147483647\n");
    std::cerr << "[resume] Direct BLAST seek to OID " << state.oid << "; validating " << state.oid_rows
              << " saved defline rows at that OID, then continuing\n";
    return prefix.string();
}

Record parse_blast_export(const std::string& line, const Taxonomy& taxonomy, bool with_oid) {
    if (!with_oid) return parse_record(line, taxonomy);
    const auto tab = line.find('\t');
    if (tab == line.npos) fail("Malformed BLAST row: expected OID followed by sequence fields");
    auto record = parse_record(std::string_view(line).substr(tab + 1), taxonomy);
    record.oid = number(std::string_view(line).substr(0, tab), "BLAST OID");
    record.export_line = line;
    return record;
}

// Human-readable BLAST metadata is only an optional progress hint. Formats we
// cannot recognize must not prevent a valid scan. BLAST may report approximate
// nucleotide lengths; a count exceeding the hint disables the estimate.
std::optional<std::uint64_t> blast_total_bases(const std::string& metadata) {
    std::istringstream lines(metadata);
    std::optional<std::uint64_t> result;
    std::string line;
    while (std::getline(lines, line)) {
        const auto fields = words(line);
        if (fields.size() != 5 || (fields[1] != "sequences;" && fields[1] != "sequence;") ||
            fields[3] != "total" || fields[4] != "bases") continue;
        auto ungroup = [](std::string_view text) {
            std::string digits(text);
            digits.erase(std::remove(digits.begin(), digits.end(), ','), digits.end());
            return number(digits, "BLAST metadata count");
        };
        try {
            if (!ungroup(fields[0])) return std::nullopt;
            const auto bases = ungroup(fields[2]);
            if (!bases || result) return std::nullopt;
            result = bases;
        } catch (const std::runtime_error&) { return std::nullopt; }
    }
    return result;
}

// ---------- 8. The complete workflow ----------

void run(Options options) {
    // Refuse even an existing empty directory: this prevents old SUCCESS files
    // or stale markers from ever being mistaken for the current run's results.
    if (options.resume) {
        if (!fs::is_directory(options.output)) fail("Resume directory does not exist: " + options.output.string());
        if (fs::exists(options.output / "SUCCESS")) fail("Run is already complete; refusing to resume or overwrite results");
    } else if (!fs::create_directory(options.output)) fail("Output directory already exists: " + options.output.string());
    RunLock lock(options.output);
    if (options.resume && fs::exists(options.output / "SUCCESS")) fail("Run is already complete");
    Manifest manifest(options.output);
    RunState state;
    Candidates candidates;
    candidates.max_load_factor(0.8F);
    std::vector<Source> sources;
    std::unique_ptr<Scanner> scanner;
    if (options.resume) {
        read_checkpoint(options, state, candidates, sources, manifest);
        std::cerr << "[resume] Restored " << phase_name(state.phase) << " checkpoint after " << state.rows << " records\n";
    }
    if (options.direct_resume && options.fasta.empty())
        options.database = resolve_blast_prefix(options.database, fs::current_path(), options).string();
    const auto started = std::chrono::steady_clock::now();
    auto last_checkpoint = started;
    const auto previous_elapsed = manifest.get("elapsed_seconds").empty() ? 0 : number(manifest.get("elapsed_seconds"), "elapsed seconds");
    bool checkpoint_ready = false;
    bool replaying = options.resume && (state.rows || state.phase == 2);
    auto save_checkpoint = [&](bool force) {
        if (!options.checkpoint_seconds || !checkpoint_ready || replaying) return;
        const auto now = std::chrono::steady_clock::now();
        if (!force && now - last_checkpoint < std::chrono::seconds(options.checkpoint_seconds)) return;
        manifest.count("elapsed_seconds", previous_elapsed + std::chrono::duration_cast<std::chrono::seconds>(now - started).count());
        manifest.count("checkpoint_records", state.rows);
        manifest.set("checkpoint_phase", phase_name(state.phase));
        write_checkpoint(options, state, candidates, sources, scanner.get(), manifest);
        manifest.save();
        last_checkpoint = std::chrono::steady_clock::now();
    };
    manifest.set("status", "running");
    manifest.set("clademers_version", VERSION);
    manifest.set("mode", "extraction");
    manifest.set("results_format", "1");
    manifest.set("lca_saved", "yes");
    if (!options.resume) manifest.set("started_utc", utc_now());
    manifest.set("session_started_utc", utc_now());
    const auto resumes = manifest.get("resume_count").empty() ? 0 : number(manifest.get("resume_count"), "resume count");
    manifest.count("resume_count", resumes + (options.resume ? 1 : 0));
    manifest.count("checkpoint_seconds", options.checkpoint_seconds);
    manifest.set("resume_mode", options.direct_resume ? "direct" : "replay");
    manifest.set("report_descendants", options.report_descendants ? "yes" : "no");
    manifest.set("write_separate_fastas", options.write_separate_fastas ? "yes" : "no");
    manifest.set("separate_fasta_compression", options.gzip_separate_fastas ? "gzip" : "none");
    manifest.set("descendant_assignment", "LCA of all observed principal and leaf taxids; one most-specific clade per retained k-mer");
    if (options.direct_resume) manifest.set("sequence_snapshot_check", "canonical file paths, sizes, mtime_ns and ctime_ns; immutable references required");
    manifest.set("zlib_version", zlibVersion());
    manifest.set("sequence_source", options.fasta.empty() ? "blast" : "fasta");
    if (options.fasta.empty()) manifest.set("database", options.database);
    else {
        manifest.set("fasta_file", fs::absolute(options.fasta).string());
        manifest.set("fasta_taxid_source", options.taxid_map.empty() ? "headers" : "map");
    }
    manifest.set("taxonomy_source", options.blast_tax.empty() ? "dmp" : "blast");
    if (options.blast_tax.empty()) {
        manifest.set("taxnodes", fs::absolute(options.nodes).string());
        manifest.set("taxnames", fs::absolute(options.names).string());
    } else {
        manifest.set("blasttax", fs::absolute(options.blast_tax).string());
        manifest.set("sqlite3", options.sqlite3);
    }
    manifest.count("target_taxid", options.target);
    manifest.count("k", options.k);
    manifest.count("scan_workers", options.threads);
    manifest.count("progress_seconds", options.progress_seconds);
    manifest.set("progress_scope", "named operation; ETA is approximate, not whole-run completion");
    manifest.set("specificity", "Exact canonical ACGT k-mers; unknown/mixed taxonomy is background; supplied reference snapshot only");
    manifest.save();
    try {
        std::cerr << "[1/4] Reading taxonomy\n";
        Taxonomy taxonomy(options);
        for (const auto& item : candidates)
            if (!taxonomy.clades.count(item.second.specificity_taxid))
                fail("Invalid checkpoint candidate specificity taxid");
        manifest.set("target_name", taxonomy.name);
        manifest.set("target_name_source", taxonomy.name_source);
        manifest.count("target_taxids", taxonomy.targets.size());
        manifest.count("taxonomy_nodes", taxonomy.node_count);
        if (!taxonomy.sqlite_version.empty()) manifest.set("sqlite3_version", taxonomy.sqlite_version);
        for (const auto& input : taxonomy.inputs) {
            manifest.set(input.key + "_file", input.path.string());
            manifest.set(input.key + "_fnv1a64", input.fingerprint.hex());
            manifest.count(input.key + "_bytes", input.fingerprint.bytes);
            manifest.set(input.key + "_fingerprint_scope", input.text ? "decompressed_content" : "file_bytes");
            if (input.text) manifest.set(input.key + "_compression", input.compressed ? "gzip" : "none");
            if (input.compressed) manifest.set("gzip", options.gzip);
        }
        std::vector<Taxid> taxids(taxonomy.targets.begin(), taxonomy.targets.end());
        std::sort(taxids.begin(), taxids.end());
        std::string taxid_text;
        for (Taxid id : taxids) taxid_text += std::to_string(id) + '\n';
        write_text(options.output / "target_taxids.txt", taxid_text);
        const auto snapshot = write_taxonomy_snapshot(options.output, taxonomy);
        manifest.set("target_taxonomy_fnv1a64", snapshot.hex());
        manifest.count("target_taxonomy_bytes", snapshot.bytes);
        std::cerr << "      " << taxonomy.name << " (" << options.target << "), "
                  << taxids.size() << " taxids including descendants\n";

        std::vector<std::string> info_args;
        std::string database_info;
        std::unique_ptr<FastaTaxidMap> mapping;
        if (options.fasta.empty()) {
            manifest.set("blastdbcmd_version", capture({options.blastdbcmd, "-version"}, options.output));
            const auto blast_help = capture({options.blastdbcmd, "-help"}, options.output);
            if (blast_help.find("-no_taxid_expansion") == std::string::npos || blast_help.find("%X") == std::string::npos)
                fail("blastdbcmd needs -no_taxid_expansion and %X support; use BLAST+ 2.15 or newer");
            info_args = blast_arguments(options); info_args.push_back("-info");
            database_info = capture(info_args, options.output);
            write_text(options.output / "database_info.txt", database_info);
        } else if (!options.taxid_map.empty()) {
            mapping = std::make_unique<FastaTaxidMap>(options);
            manifest.set("taxidmap_file", mapping->path.string());
            manifest.set("taxidmap_compression", mapping->compressed ? "gzip" : "none");
            manifest.set("taxidmap_fingerprint_scope", "decompressed_content");
            if (mapping->compressed) manifest.set("gzip", options.gzip);
            manifest.set("taxidmap_fnv1a64", mapping->fingerprint.hex());
            manifest.count("taxidmap_bytes", mapping->fingerprint.bytes);
            manifest.count("taxidmap_ids", mapping->entries.size());
        }
        manifest.save();

        const auto identity = input_identity(options, taxonomy, mapping.get(), database_info, manifest.get("blastdbcmd_version"));
        if (options.resume && state.identity != identity)
            fail("Checkpoint inputs changed: reference file paths/sizes/timestamps, taxonomy, mapping, compression, or BLAST metadata/version");
        state.identity = identity;
        checkpoint_ready = true;
        if (!options.resume) save_checkpoint(true);

        // Replay mode validates the consumed prefix. Direct BLAST validates only
        // its saved OID group; direct FASTA restores a file cursor and lookahead.
        std::uint64_t seen_rows = 0, resume_rows = state.rows;
        Fingerprint resume_prefix = state.prefix;
        std::unique_ptr<Progress> replay_progress;
        auto start_replay = [&] {
            seen_rows = 0; resume_rows = state.rows; resume_prefix = state.prefix;
            if (options.direct_resume && !options.fasta.empty()) {
                seen_rows = state.rows; replaying = false;
                return;
            }
            if (options.direct_resume && resume_rows) {
                if (!state.oid_rows || state.oid_rows > state.rows) fail("Missing saved BLAST OID group");
                seen_rows = state.rows - state.oid_rows; resume_prefix = state.oid_prefix;
            }
            if (resume_rows) {
                replay_progress = std::make_unique<Progress>(options.direct_resume ? "Validate saved BLAST OID group" :
                    "Replay saved input prefix (no k-mer calculations)", "records", resume_rows, options.progress_seconds);
                if (options.direct_resume) replay_progress->resume_at(seen_rows);
            }
        };
        auto skip_saved = [&](const Fingerprint& prefix) {
            ++seen_rows;
            if (seen_rows > resume_rows) return false;
            replay_progress->update(seen_rows);
            if (seen_rows == resume_rows) {
                if (prefix.value != resume_prefix.value || prefix.bytes != resume_prefix.bytes)
                    fail("Input prefix differs from checkpoint; results are not accepted");
                replaying = false;
                replay_progress->finish(seen_rows);
                std::cerr << (options.direct_resume ? "[resume] Saved OID group verified; continuing new work\n"
                                                   : "[resume] Saved prefix verified; continuing new work\n");
            }
            check_stop();
            return true;
        };

        std::cerr << "[2/4] Loading distinct target k-mers\n";
        auto& target_rows = state.target_rows;
        auto& target_bases = state.target_bases;
        auto commit_record = [&](const Record& record, const Fingerprint& prefix, const FastaCursor* cursor) {
            state.rows = seen_rows; state.prefix = prefix;
            if (options.direct_resume && cursor) state.fasta_cursor = *cursor;
            if (options.direct_resume && options.fasta.empty()) {
                if (record.oid >= INT32_MAX || (state.oid_rows && record.oid < state.oid))
                    fail("Direct BLAST resume requires ordered OIDs in 0..2147483646");
                if (!state.oid_rows || record.oid != state.oid) {
                    state.oid = record.oid; state.oid_rows = 0; state.oid_prefix = Fingerprint{};
                    state.blast_bases += record.sequence.size(); state.oid_length = record.sequence.size();
                } else if (state.oid_length != record.sequence.size()) fail("Inconsistent sequence length for saved BLAST OID");
                ++state.oid_rows;
                state.oid_prefix.add(record.export_line); state.oid_prefix.add("\n");
            }
            save_checkpoint(false);
            check_stop();
        };
        auto collect_target = [&](const Record& record, const Fingerprint& prefix, const FastaCursor* cursor) {
            if (skip_saved(prefix)) return;
            if (record.any_target) {
                ++target_rows; target_bases += record.sequence.size();
                const std::size_t source_index = sources.size();
                sources.push_back({record.accession, record.taxids, record.leaf_taxids, record.input_record});
                each_kmer(record.sequence, options.k, [&](const Kmer& key, std::size_t position, bool reverse) {
                    // The first target occurrence supplies representative provenance.
                    const auto inserted = candidates.try_emplace(key,
                        Origin{source_index, position, candidates.size(), reverse, record.specificity_taxid});
                    // Every observation contributes, even when the key was already
                    // present. This is independent of representative provenance.
                    if (!inserted.second)
                        inserted.first->second.specificity_taxid = taxonomy.lca(
                            inserted.first->second.specificity_taxid, record.specificity_taxid);
                });
            }
            commit_record(record, prefix, cursor);
        };
        auto& first_scan = state.first_scan;
        if (!state.phase) {
            start_replay();
            if (!options.fasta.empty()) {
                first_scan = read_fasta(options, taxonomy, mapping.get(),
                    [&](const Record& record, const Fingerprint& prefix, const FastaCursor& cursor) { collect_target(record, prefix, &cursor); },
                    "FASTA target pass", std::nullopt,
                    options.direct_resume && state.rows ? &state.fasta_cursor : nullptr, &state.prefix);
                manifest.set("fasta_compression", first_scan.compressed ? "gzip" : "none");
                manifest.set("fasta_reader", first_scan.compressed ? (options.direct_resume ? "indexed_zlib" : "gzip_command") : "plain_file");
                manifest.set("fasta_fingerprint_scope", "decompressed_content");
                if (first_scan.compressed && !options.direct_resume) manifest.set("gzip", options.gzip);
                manifest.set("fasta_fnv1a64", first_scan.fingerprint.hex());
                manifest.count("fasta_bytes", first_scan.fingerprint.bytes);
                manifest.count("fasta_records", first_scan.records);
                manifest.count("fasta_unassigned_taxid_records", first_scan.unknown_records);
                manifest.save();
                std::cerr << "      FASTA records: " << first_scan.records << "; with missing/0/N/A taxids: "
                          << first_scan.unknown_records << '\n';
            } else {
                auto target_args = blast_arguments(options, blast_resume_alias(options, state));
                target_args.insert(target_args.end(), {"-taxidlist", (options.output / "target_taxids.txt").string(),
                    "-no_taxid_expansion", "-outfmt", options.direct_resume ? "%o\t%a\t%T\t%X\t%s" : "%a\t%T\t%X\t%s"});
                Progress progress("BLAST target loading", "target rows", std::nullopt, options.progress_seconds);
                Process target(target_args, options.output);
                std::string line;
                Fingerprint prefix = options.direct_resume ? state.prefix : Fingerprint{}, anchor;
                const auto anchor_oid = state.oid;
                while (target.getline(line)) {
                    auto record = parse_blast_export(line, taxonomy, options.direct_resume);
                    const bool validating = options.direct_resume && seen_rows < resume_rows;
                    if (validating && record.oid != anchor_oid) fail("Saved BLAST OID group missing from resume export");
                    auto& current = validating ? anchor : prefix;
                    current.add(line); current.add("\n");
                    collect_target(record, current, nullptr);
                    progress.update(target_rows);
                }
                target.finish();
                progress.finish(target_rows);
            }
            if (seen_rows < resume_rows) fail("Input ended before checkpoint position");
            if (candidates.empty()) fail("No target ACGT k-mers: check target taxonomy, k, and reference coverage");
            manifest.count("target_defline_rows", target_rows);
            manifest.count("target_bases", target_bases);
            manifest.count("candidate_kmers", candidates.size());
            manifest.save();
            state.phase = 1; state.rows = 0; state.prefix = Fingerprint{};
            state.oid = UINT64_MAX; state.oid_rows = state.blast_bases = state.oid_length = 0;
            state.oid_prefix = Fingerprint{}; state.fasta_cursor = FastaCursor{};
            save_checkpoint(true);
        } else std::cerr << "      Using " << candidates.size() << " saved candidates; target loading already complete\n";
        std::cerr << "      " << candidates.size() << " candidates; full sequence identities stored in RAM\n";

        std::cerr << "[3/4] Scanning the complete reference input with " << options.threads << " workers\n";
        auto& all_rows = state.rows;
        auto& background_rows = state.background_rows;
        auto& background_bases = state.background_bases;
        auto& background_windows = state.background_windows;
        scanner = std::make_unique<Scanner>(candidates, options.k, options.threads, state.flags);
        state.flags.clear(); state.flags.shrink_to_fit();
        auto scan_background = [&](const Record& record, const Fingerprint& prefix, const FastaCursor* cursor) {
            if (skip_saved(prefix)) return;
            if (!record.only_target) {
                ++background_rows; background_bases += record.sequence.size();
                if (record.sequence.size() >= options.k) ++background_windows;
                scanner->add(record.sequence);
            }
            commit_record(record, prefix, cursor);
        };
        if (state.phase == 1) {
            start_replay();
            if (!options.fasta.empty()) {
                const auto second_scan = read_fasta(options, taxonomy, mapping.get(),
                    [&](const Record& record, const Fingerprint& prefix, const FastaCursor& cursor) { scan_background(record, prefix, &cursor); },
                    "FASTA background input", first_scan.fingerprint.bytes,
                    options.direct_resume && state.rows ? &state.fasta_cursor : nullptr, &state.prefix);
                if (first_scan.fingerprint.value != second_scan.fingerprint.value ||
                    first_scan.fingerprint.bytes != second_scan.fingerprint.bytes ||
                    first_scan.compressed != second_scan.compressed)
                    fail("FASTA input changed between passes; results are not accepted");
                manifest.set("fasta_passes_match", "yes");
                std::cerr << "      Finishing background workers\n";
                scanner->finish();
            } else {
                auto all_args = blast_arguments(options, blast_resume_alias(options, state));
                // Count each OID's sequence once for progress, even if it has multiple
                // exported deflines. Every row still goes through taxonomy/subtraction.
                all_args.insert(all_args.end(), {"-entry", "all", "-outfmt", "%o\t%a\t%T\t%X\t%s"});
                const auto total = blast_total_bases(database_info);
                if (total) manifest.count("blast_progress_total_bases", *total);
                Progress progress("BLAST background scan", "sequence bases", total, options.progress_seconds);
                if (options.direct_resume && state.rows) progress.resume_at(state.blast_bases);
                Process background(all_args, options.output);
                // Initialize the value even before the first row; this avoids
                // GCC's maybe-uninitialized diagnostic for optional's payload.
                bool have_last_oid = options.direct_resume && state.rows != 0;
                std::uint64_t last_oid = have_last_oid ? state.oid : 0;
                std::uint64_t scanned_bases = options.direct_resume ? state.blast_bases : 0;
                std::size_t last_length = options.direct_resume ? state.oid_length : 0;
                std::string line;
                Fingerprint prefix = options.direct_resume ? state.prefix : Fingerprint{}, anchor;
                const auto anchor_oid = state.oid;
                while (background.getline(line)) {
                    const auto record = parse_blast_export(line, taxonomy, true);
                    const auto oid = record.oid;
                    const bool validating = options.direct_resume && seen_rows < resume_rows;
                    if (validating && oid != anchor_oid) fail("Saved BLAST OID group missing from resume export");
                    auto& current = validating ? anchor : prefix;
                    current.add(line); current.add("\n");
                    scan_background(record, current, nullptr);
                    if (!have_last_oid || oid != last_oid) {
                        // Ordered OIDs allow constant-memory counting, also for aliases
                        // with gaps. Unusual export ordering simply removes the ETA.
                        if (have_last_oid && oid < last_oid) progress.unknown_total();
                        scanned_bases += record.sequence.size();
                    } else if (record.sequence.size() != last_length) progress.unknown_total();
                    last_oid = oid;
                    have_last_oid = true;
                    last_length = record.sequence.size();
                    progress.update(scanned_bases);
                }
                background.finish(); // A failed producer invalidates even an apparently complete stream.
                std::cerr << "      Finishing background workers\n";
                scanner->finish();
                progress.finish(scanned_bases);
            }
            if (seen_rows < resume_rows) fail("Input ended before checkpoint position");
        } else {
            std::cerr << "      Background scan already complete; rechecking inputs before output\n";
            scanner->finish();
            if (!options.fasta.empty() && !options.direct_resume) {
                const auto current = fingerprint_text(options.fasta, options, first_scan.fingerprint.bytes);
                if (current.compressed != first_scan.compressed || current.fingerprint.bytes != first_scan.fingerprint.bytes ||
                    current.fingerprint.value != first_scan.fingerprint.value)
                    fail("FASTA input changed since checkpoint");
            }
        }
        if (!background_rows || !background_windows)
            fail("No background records of length >= k; refusing an uninformative specificity result");
        if (options.fasta.empty()) {
            const auto final_info = capture(info_args, options.output);
            if (database_info != final_info) fail("Database metadata changed during scan; results are not accepted");
        }
        if (mapping) {
            mapping->verify_unchanged(options);
            manifest.set("taxidmap_passes_match", "yes");
        }
        taxonomy.verify_unchanged(options);
        if (options.direct_resume && state.identity != input_identity(options, taxonomy, mapping.get(), database_info, manifest.get("blastdbcmd_version")))
            fail("Reference files changed during run; direct-resume snapshot is no longer valid");
        manifest.set("taxonomy_inputs_match", "yes");
        manifest.count("all_defline_rows", all_rows);
        manifest.count("background_defline_rows", background_rows);
        manifest.count("background_bases", background_bases);
        replaying = false;
        state.phase = 2;
        save_checkpoint(true);
        check_stop();

        std::cerr << "[4/4] Sorting and writing retained markers\n";
        // Pointer sorting adds one pointer per retained k-mer, not a second map.
        // DNA lexical order makes output stable across thread counts/hash layouts.
        std::vector<const Candidates::value_type*> retained;
        Progress selecting("Collect retained markers", "candidates", candidates.size(), options.progress_seconds);
        std::uint64_t examined = 0;
        for (const auto& item : candidates) {
            if (!scanner->excluded(item.second.index)) retained.push_back(&item);
            // Batch the clock check for this otherwise very cheap loop.
            if (++examined % 65536 == 0) { selecting.update(examined); check_stop(); }
        }
        selecting.finish(examined);
        Progress sorting("Sort retained markers", "markers", std::nullopt, options.progress_seconds);
        std::sort(retained.begin(), retained.end(), [](const auto* a, const auto* b) { return a->first < b->first; });
        sorting.finish(retained.size());
        check_stop();
        Progress writing("Write retained markers", "markers", retained.size(), options.progress_seconds);
        std::ofstream fasta(options.output / "markers.partial.fasta");
        std::ofstream table(options.output / "markers.partial.tsv");
        if (!fasta || !table) fail("Cannot open marker outputs");
        std::unique_ptr<TaxonFastaWriter> separate;
        if (options.write_separate_fastas)
            separate = std::make_unique<TaxonFastaWriter>(options.output, taxonomy, options.gzip_separate_fastas);
        table << "marker_id\tsequence\tsource_accession\tsource_taxids\tsource_leaf_taxids\tstart_1based\tstrand";
        if (!options.fasta.empty()) table << "\tsource_record_1based";
        table << "\tspecificity_taxid\tspecificity_name";
        table << '\n';
        TaxonCounts assigned;
        Fingerprint fasta_fingerprint;
        std::uint64_t marker_index = 0;
        for (const auto* item : retained) {
            const auto& origin = item->second;
            const auto& source = sources.at(origin.source);
            const auto sequence = decode(item->first, options.k);
            const auto id = "km_" + std::to_string(options.target) + "_k" + std::to_string(options.k)
                            + "_" + std::to_string(++marker_index);
            const auto fasta_record = '>' + id + " taxid=" + std::to_string(origin.specificity_taxid) + '\n' + sequence + '\n';
            fasta << fasta_record;
            fasta_fingerprint.add(fasta_record);
            if (separate) separate->write(origin.specificity_taxid, fasta_record);
            table << id << '\t' << sequence << '\t' << source.accession << '\t' << source.taxids
                  << '\t' << source.leaf_taxids << '\t' << origin.position + 1
                  << '\t' << (origin.reverse ? '-' : '+');
            if (!options.fasta.empty()) table << '\t' << source.input_record;
            table << '\t' << origin.specificity_taxid << '\t' << taxonomy.clades.at(origin.specificity_taxid).name;
            if (options.report_descendants) ++assigned[origin.specificity_taxid];
            table << '\n';
            if (marker_index % 1024 == 0) { writing.update(marker_index); check_stop(); }
        }
        fasta.close(); table.close();
        if (!fasta || !table) fail("Marker output write failed (disk full?)");
        if (options.report_descendants)
            write_descendant_summary(options.output, options.target, taxonomy, assigned);
        if (separate) {
            separate->finish();
            manifest.count("separate_fasta_files", separate->file_count());
            manifest.set("separate_fasta_assignment", "one file per exact LCA; ancestors do not duplicate markers");
        }
        manifest.set("markers_fasta_fnv1a64", fasta_fingerprint.hex());
        manifest.count("markers_fasta_bytes", fasta_fingerprint.bytes);
        check_stop();
        fs::rename(options.output / "markers.partial.fasta", options.output / "markers.fasta");
        fs::rename(options.output / "markers.partial.tsv", options.output / "markers.tsv");
        writing.finish(marker_index);
        manifest.count("removed_kmers", candidates.size() - retained.size());
        manifest.count("retained_kmers", retained.size());
        manifest.count("elapsed_seconds", previous_elapsed + std::chrono::duration_cast<std::chrono::seconds>(
            std::chrono::steady_clock::now() - started).count());
        manifest.set("finished_utc", utc_now());
        manifest.set("status", "complete");
        manifest.save();
        write_text(options.output / "SUCCESS.partial", "Reference input and exact k-mer scan both completed.\n");
        fs::rename(options.output / "SUCCESS.partial", options.output / "SUCCESS");
        std::cerr << "Retained " << retained.size() << " of " << candidates.size()
                  << " distinct canonical k-mers. Results: " << options.output << '\n';
    } catch (const PauseRequested&) {
        save_checkpoint(true);
        manifest.set("status", "paused");
        manifest.set("paused_utc", utc_now());
        manifest.save();
        if (fs::is_regular_file(options.output / "checkpoint" / "state.bin"))
            std::cerr << "Paused. Resume with: clademers -resume " << std::quoted(options.output.string()) << '\n';
        else std::cerr << "Paused before a checkpoint was available. Restart with a new output directory.\n";
        throw;
    } catch (const std::exception& e) {
        manifest.set("status", "failed");
        manifest.set("error", e.what());
        manifest.set("finished_utc", utc_now());
        try { manifest.save(); } catch (...) { /* Preserve the original error, e.g. disk full. */ }
        throw;
    }
}

int main(int argc, char** argv) {
    try {
        if (argc == 1 || (argc == 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help"))) {
            help(); return 0;
        }
        if (argc == 2 && std::string(argv[1]) == "--version") {
            std::cout << "CladeMers " << VERSION << '\n'; return 0;
        }
        struct sigaction action{};
        action.sa_handler = request_stop;
        sigemptyset(&action.sa_mask);
        for (const int signal : {SIGUSR1, SIGTERM, SIGINT})
            if (sigaction(signal, &action, nullptr) < 0) fail("Cannot install checkpoint signal handler");
        const auto options = parse_options(argc, argv);
        if (!options.results.empty()) report_saved_results(options);
        else run(options);
        return 0;
    } catch (const PauseRequested&) {
        return 75;
    } catch (const std::bad_alloc&) {
        std::cerr << "ERROR: insufficient RAM. Target k-mers and any FASTA taxid map must fit in memory.\n";
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << '\n';
    }
    return 1;
}
