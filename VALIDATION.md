# Validation of CladeMers

## GitHub repository preparation (November 22, 2026)

The repository preparation retains version 0.13.0 and its behavior. The source
change adds only the MIT copyright/SPDX header; Makefile, tests and synthetic
data remain unchanged. Documentation now includes a quick start, portable file
links, contribution instructions, citation metadata, and publication steps.

- A clean copy containing only the 40 intended repository files compiled without
  warnings on macOS ARM64. `make test-blast` with BLAST+ 2.17.0 ran **92 tests:
  88 passed and four optional cases skipped** (one real NCBI taxonomy snapshot,
  one 0.11.0 executable fixture, and two second-architecture fixture cases).
  Earlier 0.13.0 validation with the additional executable fixtures is recorded
  below; no new algorithm or format change was introduced.
- The README quick start reproduced all six bundled synthetic markers, the
  descendant summary, and five compressed taxon FASTAs with identical decoded
  records. Python is needed for tests, not for this example or the program.
- GitHub Actions YAML passed actionlint 1.7.12; embedded shell steps and Slurm
  examples passed Bash syntax checks. Citation YAML parsed successfully. Relative
  file links resolve, and public files contain no workstation-specific paths.
- The workflow targets native Linux AMD64 and ARM64 with GCC and BLAST+ 2.17.0.
  These hosted jobs have **not yet run**; their configuration is not evidence of
  a successful Linux build. Their first execution follows the initial push.
- The initial file list excludes local research drafts, session notes, database
  data, executables, caches and historical ZIPs. Earlier archive checksums remain
  unchanged.

## Version 0.13.0: compressed separate FASTAs (November 22, 2026)

- Full optimized ARM64 suite: **92 tests, 91 passed and one skipped**. The
  optional downloaded NCBI taxonomy snapshot was unavailable. Actual BLAST+
  2.17.0 and 0.11.0 saved-result compatibility tests passed. ARM64 and x86-64
  builds compiled without warnings.
- All **25 descendant/export tests passed** under AddressSanitizer plus
  UndefinedBehaviorSanitizer, and on x86-64 through Rosetta. No sanitizer
  diagnostics occurred; unsupported macOS leak detection was disabled.
- Four new tests cover `-wsf-gz` / `--write-separate-fastas-gz`, implied `-wsf`,
  either flag order, unchanged combined catalogs, independent `-rd`, four-file
  portable export with original inputs deleted and PATH empty, checkpoint
  persistence, and late enablement of compression during output recovery.
- Both plain and gzip outputs pass empty-catalog, corrupt-catalog rejection,
  actual 0.11.0 result export, and 129-taxon stress cases under a 64-descriptor
  limit. Decoding reproduces the original FASTA records and within-file order
  exactly. The stress fixture verifies multiple gzip members after cache
  eviction/reopening, all CRCs/trailers, and `gzip -t` acceptance.
- Target/background SIGKILL recovery restores compression without repeating
  the flag. These checkpoints were exchanged between ARM64 and x86-64 builds
  in both directions. Output recovery replaces stale plain and compressed
  generated files without mixing formats or duplicating records.
- A file-size limit forces gzip finalization to fail while writing its trailer.
  The error is reported without publishing final FASTAs, index or SUCCESS.

Checkpoint state format is 5 to save the compression setting; candidate format
remains 2 with 56-byte entries. Resume requires the same 0.13.0 release. Portable
results remain format 1, compatible with completed 0.11.0 and newer runs.
Compression uses existing zlib at level 6, without another dependency or an
external gzip process. No GNU GCC/Linux, full NT throughput, large HPC export
benchmark, or native AMD/Intel hardware test was performed in this update.

## Version 0.12.0: separate FASTAs by LCA (November 22, 2026)

- Full optimized ARM64 suite: **88 tests, 87 passed and one skipped**. The
  optional downloaded NCBI taxonomy snapshot was unavailable. Tests included
  real BLAST+ 2.17.0 databases, SQLite/gzip input, and an actual 0.11.0 executable
  to generate prior-version results for export. Default compilation was clean.
- All **21 descendant/export tests passed** under AddressSanitizer plus
  UndefinedBehaviorSanitizer, and on x86-64 through Rosetta. No sanitizer
  diagnostics occurred; unsupported macOS leak detection was disabled.
- Seven new tests verify the two aliases, independence from `-rd`, exact
  partitions of the combined marker catalog, preserved record bytes/order,
  numeric filenames, indexed names and counts, root assignments, no ancestor
  duplication, zero-count taxa, and empty catalogs.
- A stress fixture produces 129 nonempty taxon files while limiting the process
  to 64 file descriptors. Extraction and saved-result export both pass, checking
  cache eviction/reopening without truncation or duplicate appends. Scientific
  names include duplicates, punctuation, a slash, and Unicode; filenames remain
  numeric. An independent expected sequence-to-taxid mapping checks assignments.
- Saved-result export works with four moved files, original inputs deleted and
  PATH empty. Real 0.11.0 results work without conversion or another extraction.
  A changed saved catalog fails fingerprint validation without publishing final
  FASTAs, index, summary or SUCCESS.
- Target/background SIGKILL recovery restores `-wsf` without repeating it.
  Checkpoints with the flag enabled were exchanged in both directions between
  ARM64 and x86-64 builds. Enabling `-wsf` at output-phase resume replaces stale
  generated files rather than appending duplicates. Unexpected files in a
  staging directory are preserved and cause a clear failure.
- The bundled preview now includes five separate FASTAs containing six markers,
  plus `taxon_fastas.tsv`. The preview test verifies exact output bytes.

Checkpoint state format is 4 to store the separate-FASTA option. Candidate
format remains 2 with 56-byte entries; 0.12.0 requires its own checkpoint
version. Portable results format remains 1, so completed 0.11.0 extractions
remain compatible. GNU GCC/Linux, full NT throughput, HPC filesystem overhead
for large exports, and native AMD/Intel hardware were not tested here.

## Version 0.11.0 report column naming update (November 22, 2026)

The summary column `clade_specific_kmers` is now named
`assigned_kmers_incl_descendants`. It still counts assignments to the node plus
all descendants; `assigned_kmers` still counts assignments directly to the node.
The C++ change is limited to the output header string. Version, calculations,
checkpoint serialization and portable marker catalogs remain 0.11.0-compatible.

The default build completed without warnings. Three existing focused tests
passed: independent assignments versus separate descendant runs, deferred
reporting without original inputs, and exact bundled-preview output. The source
package is `CladeMers-0.11.0-report-labels-source.zip`; the earlier archive is
preserved. Full-suite results for the underlying 0.11.0 release follow below.

## Version 0.11.0: always saved LCA and later reports (November 22, 2026)

- Full optimized ARM64 suite: **81 tests, 80 passed and one skipped**. The
  optional downloaded NCBI taxonomy snapshot was unavailable. Real BLAST+
  2.17.0 databases and synthetic binary/SQLite taxonomy fixtures were tested.
  Default compiler warning flags produced no warnings.
- All **14 descendant/report tests passed** under AddressSanitizer plus
  UndefinedBehaviorSanitizer, and on an x86-64 build through Rosetta. No
  sanitizer diagnostics occurred; unsupported macOS leak detection was disabled.
- Extraction without `-rd` saves the same LCA-labeled marker FASTA and TSV as
  extraction with `-rd`. The independent all-occurrence string/set oracle,
  separate per-descendant runs, bit boundaries, threading, unknown/mixed taxids,
  broad annotations, both taxonomy inputs, and compressed text inputs still pass.
- Deferred reports match immediate reports byte for byte. A test copies only
  `markers.fasta`, `target_taxonomy.tsv`, `manifest.tsv`, and `SUCCESS` to another
  directory, deletes the original database/taxonomy/reader, and empties PATH.
  Reporting still succeeds without marker TSV or checkpoint files and leaves
  all four source files unchanged. A separate test removes compressed source
  FASTA and BLAST taxonomy after extraction and generates the expected report.
- Tests reject missing completion markers, incomplete and legacy manifests,
  missing/truncated/changed marker catalogs, absent or invalid LCA labels,
  changed names, disconnected/cyclic taxonomy, absent parents, duplicate taxids,
  invalid CLI combinations, and attempts to overwrite an existing report.
  Empty catalogs and root/leaf targets produce the expected zero-count reports.
- Interrupted runs now begin without `-rd`, then enable it on resume. Target,
  background, and output recovery preserve LCA labels, saved taxonomy and
  summary counts. Plain and indexed-gzip checkpoints were exchanged in both
  directions between ARM64 and x86-64 builds. Existing direct BLAST OID-range,
  replay, corruption, input-change, locking and Slurm-wrapper tests also pass.
- The bundled `examples/descendant_demo` contains synthetic input and actual
  generated preview outputs. Its test verifies six marker records and the
  deferred report: direct/subtree counts 1/6, 1/4, 2/2, 1/1, 1/1, and 0/0.

Portable results format 1 includes saved FASTA/taxonomy fingerprints and
completion metadata. Checkpoint state/candidate layouts remain 3/2 (56 bytes
per serialized candidate), but 0.11.0 requires its own checkpoint version;
older checkpoints need their original executable. Earlier outputs cannot use
the new `-results` mode. Report-only work is not checkpointed and never reloads
reference inputs. GNU GCC/Linux, full NT performance, HPC memory usage, and
native AMD/Intel hardware were not tested here.

## Version 0.10.0: descendant-specific reports (November 22, 2026)

- The full **75-test suite completed: 74 passed, one skipped**, on macOS ARM64
  with Apple Clang, BLAST+ 2.17.0, SQLite, gzip, and an x86-64 alternate build.
  The optional downloaded NCBI taxonomy snapshot was unavailable; synthetic
  SQLite/binary taxonomy fixtures and real version-5 BLAST databases were tested.
  Compilation with the default warning flags emitted no warnings.
- Eight new descendant tests compare exact marker-to-taxid assignments with an
  independent Python DNA-string/ancestor-set oracle. They cover k=1, 7, 32, 33,
  42, and 64, sibling sharing, parent-level annotations, multiple principal and
  leaf taxids, reverse complements, lowercase, mixed/unknown background, empty
  results, root/leaf targets, zero-count taxa, and deterministic threaded output.
  Inclusive summary counts also match separate runs targeting each descendant.
- FASTA header and mapping inputs, plain/indexed gzip FASTA, compressed mappings
  and dumps, and both taxonomy modes agree. Missing internal BLAST names use
  taxid labels. A real BLAST database test includes identical sequences assigned
  to different species and compares all occurrences with the same oracle.
- SIGKILL recovery in target and background phases, output regeneration, and
  plain/indexed-gzip direct resume preserve marker labels and summary bytes.
  A late record broadens a saved species assignment to its shared ancestor.
  Descendant checkpoints were exchanged in both directions between ARM64 and
  x86-64 builds; the x86-64 build ran locally through Rosetta.
- All eight new descendant tests also passed under AddressSanitizer and
  UndefinedBehaviorSanitizer, and on x86-64/Rosetta (seven in each initial suite,
  followed by the real-BLAST test once its executable was available). No
  sanitizer diagnostics occurred. Unsupported macOS leak detection was disabled.
- Without `-rd`, the marker FASTA and original TSV schema remain unchanged.
  Both aliases are tested; duplicate aliases and adding the option on resume
  are rejected. The saved option is automatically restored on resume.

Checkpoint state format is now 3 and candidate format 2; each serialized
candidate uses 56 bytes plus source-table storage. Older checkpoints require
their original CladeMers release. The new report needs a run started with `-rd`.
The 0.9.0 GCC progress-counter fix is retained. GNU GCC/Linux, full NT runs,
HPC memory/performance, and native Intel/AMD hardware were not tested here.

## Version 0.9.0 GCC warning fix (November 21, 2026)

The user reported `-Wmaybe-uninitialized` for the BLAST progress counter's
`std::optional<uint64_t>` payload when compiling with g++ on a Linux AMD64 server.
Every dereference was guarded by the optional's presence check. The source now
uses an explicitly initialized integer and a separate presence flag. OID zero,
repeated OIDs, ordering checks, and the saved OID on resume retain their behavior.
No warning flags are disabled. Program version and checkpoint format remain
0.9.0; only progress-counter state representation changed.

The patched source builds without warnings with Apple Clang 21 on macOS ARM64.
Seven existing tests passed: duplicate-OID progress/ETA, missing metadata,
unreliable totals/out-of-order OIDs/failed exports, a real version-5 BLAST
database, direct resume through filtered aliases and repeated deflines,
BLAST/SQLite replay resume, and SIGKILL/repeated-pause recovery. GNU GCC is not
installed in the local environment, so the original server compiler diagnostic has
not been retested here. The updated archive is named
`CladeMers-0.9.0-gcc-fix-source.zip`; the original 0.9.0 archive is preserved.

## Version 0.9.0: direct resume and indexed gzip (November 21, 2026)

- All **67 tests passed** with the optimized native ARM64 build: the full
  66-test suite plus the subsequently added concatenated-member seek test.
  ASan+UBSan passed 65 tests with the two optional architecture-exchange tests
  skipped. Builds emitted no warnings and tests produced no sanitizer diagnostics.
  Leak detection is unsupported by this macOS sanitizer and was disabled.
- All eight direct-resume tests passed with the macOS x86-64 build under Rosetta.
  Indexed gzip background checkpoints were exchanged **in both directions**
  between ARM64 and x86-64 builds. The checkpoint stores explicit integer fields
  and dictionary bytes, not a native zlib object or pointer.
- Tests interrupt real processes after periodic checkpoints and compare resumed
  marker FASTA/TSV bytes, source coordinates/record numbers, target subtree,
  counters, and full-content FASTA fingerprints with uninterrupted runs.
  Plain FASTA seeks directly; gzip resumes target and background phases from
  non-header compressed offsets. Standard-compressor fixtures require less than
  8 MiB of local decompression; this is a test observation, not a general bound.
- Gzip tests cover stored/fixed/dynamic DEFLATE blocks, optional header fields
  and header CRC, concatenated and empty members, member boundaries inside
  sequence lines, CRLF, missing final newline, padding, corruption and truncation.
  Direct resume after multiple gzip members also matches uninterrupted results.
- A real version-5 BLAST database split across two database prefixes and nested,
  range-filtered aliases resumes both phases. An instrumented real reader repeats
  deflines to exercise a partially consumed OID group. Its audit confirms that
  resume begins at the saved global OID and that the background export contains
  no earlier OIDs. The original alias filters, source provenance, and counters
  are preserved. Paths containing spaces and shell metacharacters remain literal.
- FASTA edits preserving size and restored mtime are detected by ctime. Alias
  and physical BLAST-volume timestamp changes independently invalidate direct
  checkpoints before sequence export. These are metadata checks for immutable
  inputs, not full checksums of a skipped reference prefix.
- Earlier algorithm, taxonomy, compressed-text, signal, corruption, locking,
  output-recovery, and local Slurm-wrapper tests still pass. Compatibility tests
  for the older external readers explicitly choose `-resume-mode replay`.

Environment: Apple Clang 21, macOS ARM64, zlib 1.2.12, Apple gzip, SQLite CLI,
and BLAST+ 2.17.0+. x86-64 executables run under Rosetta with native ARM64 reader
subprocesses. Native AMD/Intel, Linux/HPC execution, full NT performance, and
shared-filesystem restart latency remain untested. Gzip access-point spacing
depends on DEFLATE block boundaries; unusually large blocks can require longer
local replay. Taxonomy/maps and checkpoint tables still load on restart.
0.8.0 checkpoints require their original release; 0.9.0 adds a new cursor format.

## Version 0.8.0: checkpoint/resume (November 21, 2026)

- All **59 tests passed** with the optimized native ARM64 build, including 12
  restart tests. The same suite under ASan+UBSan passed 58 tests with the optional
  second-architecture test skipped; no warnings or sanitizer diagnostics.
- All 12 restart tests passed with a macOS x86-64 executable under Rosetta.
  Target-loading and background checkpoints were exchanged in **both directions**
  between native ARM64 and x86-64 builds, preserving markers, source provenance,
  and final counters. Same release, inputs and reader versions were used.
- Real SIGTERM/SIGUSR1/SIGINT pauses and process-group SIGKILL interruption were
  followed by successful resumes. Target-loading tests cover k=1/32/33/42/64.
  Saved background exclusions survive worker-count changes and records spanning
  the 1 MiB chunk boundary. A second interruption during prefix replay leaves
  the existing checkpoint unchanged.
- Gzip target/background resumes use concatenated FASTA members plus compressed
  mappings and both taxonomy dumps. A genuine version-5 BLAST database resumes
  with SQLite/binary BLAST taxonomy and a changed reader executable path returning
  the same BLAST version/output. Plain FASTA output-phase recovery is also tested.
- Resumed runs match uninterrupted marker FASTA/TSV bytes, representative
  coordinates/strand, target subtree, and candidate/background/retained counts.
  Output-phase recovery rewrites partial/final marker outputs consistently.
- Changed sequence prefixes, truncated replays, altered BLAST metadata,
  taxonomy, mapping, or completed FASTA content are rejected. Corrupt/truncated
  checkpoint state, corrupt/missing candidate tables, unsupported CLI overrides,
  completed outputs, and concurrent resume writers are rejected without SUCCESS.
- A forced checkpoint-write failure preserves the old state; recovery succeeds
  once the obstruction is removed. Abandoned `.partial` files are ignored.
  Disabled checkpoint mode produces no checkpoint directory.
- The Slurm example passes `bash -n`; a local shell integration test sends USR1,
  checks forwarding, waits for a simulated delayed save, observes exit 75, and
  verifies that the next submission selects resume. This is **not** a live Slurm
  scheduler test. Scheduler timing, shared-filesystem semantics, disk requirements,
  and full-NT replay/checkpoint performance still require HPC measurements.

Environment: Apple Clang 21 on macOS ARM64, Apple gzip, SQLite CLI, and NCBI
BLAST+ 2.17.0+. x86-64 binaries ran under Rosetta with native ARM64 reader
subprocesses. Native AMD/Intel, Linux/HPC, power-loss recovery, and real NT
performance have not been tested here. Atomic rename/fsync and advisory locking
must be supported by the result filesystem. Checkpoint checksums are FNV-1a64
accidental-damage checks, not cryptographic authentication.

## Version 0.7.0: progress percentages and ETA (November 21, 2026)

- All **47 tests passed** with the optimized macOS ARM64 build, including the
  genuine BLAST+ version-5 integration and real NCBI taxonomy snapshot tests.
- Seven focused tests also passed under AddressSanitizer + UndefinedBehaviorSanitizer
  and with a macOS x86-64 build under Rosetta: all five new progress tests, genuine
  BLAST integration, and the independent string oracle across k/worker boundaries.
  Builds reported no warnings; no sanitizer diagnostics were reported.
- Controlled, delayed BLAST exports exercise an actual intermediate percentage
  and numeric ETA. Consecutive duplicate OIDs contribute once to progress, while
  every defline still participates in subtraction. Tests cover comma-separated
  metadata counts, unknown/zero/overflowed totals, underestimated totals, OID
  ordering changes, empty retained results, and failed exports without a false
  scan-completion line. The genuine BLAST metadata format is recognized.
- Plain FASTA reports byte totals in both passes. Concatenated gzip FASTA reports
  an unknown total initially and the measured full content size on its second
  pass. Compressed maps/dumps similarly acquire totals for rechecks. Existing
  invocation-count assertions confirm no additional decompression/counting pass.
- CLI interval bounds, stderr-only output, final counts, and all previous
  fingerprint/failure/marker-equivalence tests pass. Slurm example syntax passes.
- Progress percentages describe named operations; ETAs use their average rate.
  Timing accuracy under real NT/HPC load has not been benchmarked. Unknown totals
  are explicit, and read completion does not replace final snapshot validation.

The environment remains macOS ARM64, Apple Clang 21, Apple gzip, SQLite CLI,
and NCBI BLAST+ 2.17.0+. The x86-64 executable ran under Rosetta with native
ARM64 subprocesses. Linux/HPC, native AMD/Intel, and full-NT runs remain untested.

## Version 0.6.0: streamed gzip taxonomy dumps (November 21, 2026)

- All **42 tests passed** with the optimized macOS ARM64 build, with
  AddressSanitizer + UndefinedBehaviorSanitizer, and with a macOS x86-64 build
  under Rosetta. No compiler warnings or sanitizer diagnostics were reported.
- All 16 plain/gzip combinations across FASTA, map, nodes, and names produced
  identical markers, source provenance, target labels, and descendant lists.
  Header-based FASTA without a map was also checked with compressed dumps.
- Plain/gzip nodes and names, independently and together, produced matching
  markers with a genuine version-5 BLAST database. Existing BLAST binary/SQLite
  taxonomy tests, including the real NCBI snapshot, still pass.
- Dump tests cover concatenated/empty gzip members, boundaries inside fields,
  long name records across reader buffers, CRLF, no final newline, and exact
  decompressed-content fingerprints. Corrupt/truncated headers, bodies and
  trailers, CRC failures, malformed fields, invalid/cyclic/duplicate taxonomy,
  and missing/duplicate scientific names are rejected.
- Missing readers, failed reader exits during loading/final verification,
  actual dump edits or truncation after loading, and a compression-mode change
  with identical content all invalidate the run without `SUCCESS`.
- Plain dump fingerprints retain their old meaning. The tree and name parsing
  rules and the exact k-mer algorithm are unchanged. BLAST binary/SQLite files
  retain raw-file fingerprinting and SQLite snapshot checks.
- The Makefile and compile/link dependencies are unchanged. Slurm example
  syntax was checked after adding compressed-dump guidance.

The environment is the same as below: macOS ARM64, Apple Clang 21, Apple gzip,
SQLite CLI, and NCBI BLAST+ 2.17.0+. The x86-64 executable used native ARM64
subprocesses under Rosetta. Linux/HPC execution, native AMD/Intel hardware,
GNU gzip, and full-NT resource requirements remain untested here.

## Version 0.5.0: streamed gzip taxid maps (November 21, 2026)

- All **37 tests passed** with the optimized macOS ARM64 build, with
  AddressSanitizer + UndefinedBehaviorSanitizer, and with a macOS x86-64 build
  under Rosetta. No compiler warnings or sanitizer diagnostics were reported.
- All four plain/gzip FASTA/map combinations yielded identical marker tables
  and decompressed-content fingerprints with both taxonomy modes and 1/30
  scanning workers. All previous BLAST, FASTA, and gzip tests still pass.
- New map tests cover magic-based compression detection, concatenated/empty
  gzip members, long comment lines crossing reader buffers, CRLF, missing final
  newlines, headers/comments, duplicate-ID unions, missing and unknown taxids,
  exact accession/pipe matching, and conflicting header/map declarations.
- Corrupt/truncated maps, CRC failures, invalid mapping text, unavailable
  readers, and nonzero reader exits during loading or final verification are
  rejected without `SUCCESS`. Tests also modify or truncate the compressed
  file after loading and verify that final checking rejects it. Controlled
  readers test changed taxids and final-newline status during verification.
- FASTA and maps now share a checked streaming text reader. Map verification
  rereads content without constructing a second map. The decompressed map's
  byte count/fingerprint, compression mode, and successful recheck are recorded
  in the manifest. Plain-map fingerprints retain their previous meaning.
- Compilation, Makefile, and dependencies are unchanged from 0.4.0. Gzip is
  needed when either input is compressed; neither input is decompressed to disk.
  Map entries still occupy RAM regardless of their on-disk compression.

Environment and platform limits are the same as for 0.4.0 below. These checks
do not establish Linux/HPC execution, native AMD/Intel execution, GNU gzip
behavior, or full-NT performance.

## Version 0.4.0: streamed gzip FASTA (November 21, 2026)

- All **32 tests passed** with the optimized macOS ARM64 build, with
  AddressSanitizer + UndefinedBehaviorSanitizer, and with a macOS x86-64 build
  under Rosetta. No compiler warnings or sanitizer diagnostics were reported.
- Compressed and plain FASTA produced identical markers, provenance, and
  content fingerprints using header taxids or a mapping file, both taxonomy
  modes, and 1/30 scanning workers. Earlier BLAST and FASTA tests still pass.
- Gzip tests cover content-based detection with non-`.gz` filenames,
  concatenated/empty members, member boundaries inside headers/sequences,
  CRLF, a missing final newline, long lines spanning the subprocess buffer,
  and background k-mers spanning the scanner's 1-MiB chunk boundary.
- The independent string oracle verifies compressed input for
  k=1,2,7,31,32,33,42,63,64. Malformed/truncated gzip headers, bodies and trailers,
  CRC errors, empty decoded input, and a truncated later member are rejected.
- Controlled readers verify rejection of first- or second-pass nonzero exit
  statuses even after complete-looking output, and changed content/newline
  status between passes. Missing readers fail clearly; plain FASTA does not
  execute a gzip reader. Shell metacharacters in paths remain literal.
- The Makefile and compile/link flags are unchanged. Gzip is a runtime CLI
  dependency only for compressed FASTA; it streams into the parser without
  writing a decompressed reference file. Taxid maps remain plain text.

The gzip reader used here was `/usr/bin/gzip`, reporting Apple gzip
479.0.0.700.1. The other tools and real NCBI taxonomy snapshot are the same as
below. The x86-64 CladeMers build used native ARM64 subprocesses under Rosetta.
Linux/HPC execution, native x86-64 hardware, GNU gzip, and full-NT performance
were not tested in this session. `CLADEMERS_GZIP` can select an absolute gzip
reader path when reproducing the developer tests.

## Version 0.3.0: FASTA and BLAST sequence input (November 21, 2026)

- All **27 tests passed** with the optimized macOS ARM64 build, with
  AddressSanitizer + UndefinedBehaviorSanitizer, and with a macOS x86-64 build
  under Rosetta. Builds had no compiler warnings or sanitizer diagnostics.
- FASTA with header taxids or a two-column mapping file produced identical
  marker catalogs and common provenance fields to a genuine version-5 BLAST
  database containing the same sequences and assignments. Both `.dmp` and
  BLAST taxonomy inputs were exercised, including the real NCBI snapshot below.
  FASTA intentionally adds a record-number column and uses `N/A` for the
  BLAST-specific leaf-taxid field.
- The independent DNA-string oracle covers k=1,2,7,31,32,33,42,63,64 with serial
  and parallel FASTA scans. Tests check coordinates, reverse complements,
  ambiguity/gaps, wrapped sequences, CRLF, missing final newlines, repeated IDs,
  Kraken tags, bracketed annotations, and exact raw-byte fingerprints.
- Tests cover mixed/unknown taxids, map unions, exact versioned/pipe ID matching,
  missing map entries, conflicting declarations, malformed input, empty records,
  gzip rejection, CLI exclusivity, and missing files. FASTA tests run with an
  empty executable search path: `.dmp` mode launches no external program, and
  BLAST taxonomy mode uses only an explicitly selected SQLite CLI.
- A supplemental local C++ harness synchronously changed inputs between the
  two FASTA passes. Appending a newline, a same-size FASTA edit, and a same-size
  mapping edit were each rejected without `SUCCESS`; the harness was also
  compiled with both sanitizers. It exercised `run()` through a diagnostic
  stream callback, without timing races or product test hooks.
- Existing BLAST behavior, taxonomy validation, and the exact k-mer engine
  remain covered by all 20 earlier tests. The BLAST marker TSV keeps its
  original seven columns.

Environment and reference snapshot are the same as for 0.2.0 below. The
x86-64 executable ran under Rosetta on ARM64 macOS with native ARM64 BLAST+
and SQLite subprocesses. These checks do not establish Linux execution,
native AMD/Intel execution, or NT-scale performance. FASTA currently requires
regular uncompressed files and holds the optional taxid map in memory.

## Version 0.2.0: two taxonomy input modes (November 21, 2026)

- All **20 tests passed** with the optimized macOS ARM64 build, including real
  BLAST+ database integration and an actual NCBI BLAST taxonomy snapshot.
- The same 20 tests passed with AddressSanitizer + UndefinedBehaviorSanitizer,
  without sanitizer diagnostics, and with an x86-64 build under Rosetta.
- Both `.dmp` and `-blasttax` modes produced identical marker tables for
  equivalent trees with 1 and 30 workers. A separate recursive SQL oracle
  confirmed the descendant list for the real NCBI snapshot.
- New tests cover root/leaf/internal targets, big-endian binary names,
  missing internal names with numeric-label fallback, mode exclusivity,
  missing/corrupt files, malformed SQLite schemas, invalid taxonomy trees,
  active journal/WAL detection, unavailable/failed/truncated SQLite readers,
  path quoting, provenance, and taxonomy changes before completion.
- K-mer encoding, record taxonomy classification, and parallel subtraction
  were verified unchanged from 0.1.0. Scientific names are labels only.

Environment: Apple Clang 21.0.0, Python 3.9.6 for tests, SQLite CLI 3.51.0,
and official NCBI BLAST+ 2.17.0+. The x86-64 executable was run under Rosetta
on ARM64 macOS and used the native ARM64 BLAST+ and SQLite subprocesses.

The real taxonomy test used NCBI's
[`taxdb.tar.gz`](https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz), downloaded
November 21, 2026; its three files were dated November 20, 2026. Archive SHA256:

```text
f8d1892a63b434f6afaa32f8b9e2f9acc8487e4d6db2978850a55049eea6efcc
```

Its SQLite tree contained 3,476,917 nodes and a 143-taxid Placozoa subtree.
The test derived a matching `.dmp` subtree from that snapshot and obtained the
same result in both modes: 2 synthetic candidate 42-mers, 1 retained. The
binary name index omitted taxid 10226, so the BLAST-mode label was `taxid 10226`.
A second run targeting 9606 verified the real `Homo sapiens` name record and
correctly retained zero markers because the fixture's human sequence was shared.
No NT scan or biological detection is implied by these synthetic-reference tests.

## Earlier validation: version 0.1.0

Initial validation was performed on November 20, 2026, on macOS ARM64 with
Apple Clang 21.0.0, Python 3.9.6 for the developer tests, and NCBI BLAST+ 2.17.0+
(the official ARM64 macOS binaries, build July 1, 2025).

## Results

- Compiled with `-O3 -std=c++17 -pthread -Wall -Wextra -Wpedantic`, without
  compiler warnings.
- **11 tests passed**, including a genuine version-5 nucleotide BLAST database
  constructed locally with `makeblastdb` and read using `blastdbcmd`.
- The entire 11-test suite also passed with **AddressSanitizer and
  UndefinedBehaviorSanitizer**, without sanitizer diagnostics.
- The bundled synthetic demo produced 2 candidate 42-mers and retained exactly
  1, excluding the sequence shared with background.
- The real-BLAST test produced identical marker tables using 1 and 30 workers.

The string oracle uses Python substring extraction and reverse complementation,
independently of the C++ rolling bit representation. Tested k values are
1, 2, 7, 31, 32, 33, 42, 63, and 64. Synthetic process fixtures exercise
producer failures after output, incomplete rows, mixed and unknown annotations,
taxonomy errors, changed database metadata, and paths containing shell
metacharacters. The real database is moved into a directory containing a space
to verify the BLAST reader's database-name quoting.

## Rename verification (November 21, 2026)

After renaming the project to CladeMers, the source to `clademers.cpp`, and the
executable to `clademers`:

- The optimized build completed without compiler warnings, and all **11 tests
  passed**, including the real BLAST+ version-5 integration test.
- Help/version output, the default `clademers_10226_k42` output directory, and
  the `clademers_version` manifest key were checked using the bundled demo.
  The demo again produced 2 candidate 42-mers and retained exactly 1.
- Installation into a temporary staging directory produced a working
  `bin/clademers`; the `CLADEMERS_BINARY` developer-test override and Slurm
  example syntax were also checked.
- Source and test changes were verified to be name substitutions only.
  The sanitizer results above are from November 20; sanitizers were not rerun
  for this rename.

## Reproduce

```bash
make test-blast BLAST_BIN=/path/to/ncbi-blast/bin \
  BLAST_TAXDB=/path/to/extracted/blast-taxonomy
```

`BLAST_TAXDB` is optional; omit it to skip the real NCBI taxonomy snapshot test.
The new taxonomy-mode tests need the SQLite CLI; `CLADEMERS_SQLITE3` can select
its path. Tests needing unavailable external tools/data are explicitly skipped.

To reproduce the sanitizer build with the current filenames:

```bash
c++ -O1 -g -std=c++17 -pthread -Wall -Wextra -Wpedantic \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  clademers.cpp -lz -o /tmp/clademers-asan
CLADEMERS_BINARY=/tmp/clademers-asan BLAST_BIN=/path/to/ncbi-blast/bin \
  BLAST_TAXDB=/path/to/extracted/blast-taxonomy \
  python3 tests/test_clademers.py
```

## HPC source-package verification (November 21, 2026)

- Added a phony `make clean` target for the executable, optional object file,
  and macOS debug-symbol directory. Checks in a temporary workspace confirmed
  that cleaning succeeds twice, including with a file named `clean`, removes
  only those artifacts, and preserves sources, databases, results, and archives.
- A clean native ARM64 build on macOS passed all **11 tests**, including the
  genuine BLAST+ version-5 integration test, without compiler warnings.
- The same source compiled for macOS x86-64 with Apple Clang's `-arch x86_64`
  and passed all **11 tests** under Rosetta. The BLAST+ subprocess was the
  existing native ARM64 BLAST+ 2.17.0+ installation. This is an x86-64 executable
  test under translation, not a native AMD/Intel machine or Linux test.
- The default build flags contain no CPU-specific instruction-set options.
  Native compilation is intended for AMD64/x86-64 and ARM64/AArch64 with a
  suitable C++17/POSIX toolchain. The C++ algorithm was unchanged for packaging.

## Not yet established

No full NT scan has been run with this implementation. HPC memory requirements,
throughput, thread scaling, and behavior with the user's specific NT snapshot
have not been measured. Linux execution, native x86-64 hardware execution,
and 32-bit ARM execution have not been tested in this session. No biological
marker validation or SRA search has been performed by this software.
Passing the small tests is not evidence of
NT-scale performance or organism-level specificity.
