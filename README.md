# CladeMers

**Find every distinct, exact k-mer present in a target lineage and absent from
the background in a local, taxonomy-annotated nucleotide BLAST database or FASTA.**

CladeMers is a small C++ command-line program. It reads an existing NT database
or plain/gzip-compressed FASTA with taxids in its headers or a separate mapping file.

This is version **0.13.0**. Every run saves each marker's most-specific clade
(LCA) in its FASTA header and TSV, plus a portable taxonomy snapshot.
`-rd -results DIRECTORY` generates a descendant report later without rescanning
NT or needing the original inputs. `-wsf` writes one FASTA per most-specific LCA
taxid, during extraction or later from saved results. Use `-wsf-gz` for gzip-compressed
separate FASTAs. Direct checkpoint resume
supports BLAST, plain FASTA,
and indexed gzip FASTA, with progress percentages/estimated time remaining in
the log. All text inputs can stay gzip-compressed: FASTA, taxid maps, `nodes.dmp`, and `names.dmp`. Both sequence sources work with either
NCBI `.dmp` or BLAST taxonomy inputs.
Use the original release to resume older checkpoints; portable result reporting
starts with 0.11.0. See [VALIDATION.md](VALIDATION.md) for tested cases and
limitations. Full NT throughput, memory measurements, and HPC performance
benchmarking are still needed.

## Quick start

With a C++17 compiler, Make, and zlib development headers installed, run this
small synthetic example. If you downloaded a source ZIP, unpack it and start
with `make` from that directory instead of cloning:

```bash
git clone https://github.com/deanbobo/CladeMers.git
cd CladeMers
make
./clademers -fasta examples/descendant_demo/references.fasta \
  -taxnodes examples/descendant_demo/nodes.dmp \
  -taxnames examples/descendant_demo/names.dmp \
  -taxid 10226 -k 42 -out clademers_demo

# Generate a descendant summary and compressed per-taxon FASTAs later:
./clademers -rd -wsf-gz -results clademers_demo -out clademers_demo_gzip
```

This example needs neither BLAST+ nor Python. It produces six synthetic markers;
see the [expected descendant report](examples/descendant_demo/README.md).
Use new output directory names when repeating a run.

- [HPC setup and NT/Placozoa commands](HPC.md)
- [Input and command-line options](#dependencies-and-build)
- [Developer tests](#developer-tests) and [contributing](CONTRIBUTING.md)
- [Software citation](CITATION.cff)

## License and citation

CladeMers is developed by Dean M. Bobo and distributed under the
[MIT License](LICENSE). Dependencies and reference databases retain their own
licenses and terms.

If you use CladeMers, please cite the software using [CITATION.cff](CITATION.cff)
and report the version, command, reference database release, and taxonomy
snapshot. The citation describes the software; no journal article or archived
release DOI is currently specified. Report problems in the
[GitHub issue tracker](https://github.com/deanbobo/CladeMers/issues).

## Dependencies and build

To build, use **GCC 9+ or a recent Clang with C++17 support**, on Linux or macOS,
and **zlib headers and library** (`zlib.h`, linked with `-lz`). Load your HPC's
zlib development module/package; Debian/Ubuntu calls it `zlib1g-dev`, and
RHEL/Rocky/AlmaLinux calls it `zlib-devel`. The Makefile links zlib automatically.
For **BLAST sequence input**, have **NCBI BLAST+ 2.15 or newer**, with
`blastdbcmd` on your `PATH`. The database must support taxid-list extraction
(use a version-5 database). **FASTA input does not require BLAST+.**
The **optional `-blasttax` mode also needs the `sqlite3` command-line program**
with read-only support, on `PATH` or supplied with `-sqlite3 /path/to/sqlite3`.
The `.dmp` mode does not use SQLite. Neither mode links against a SQLite library.
**Compressed maps and taxonomy dumps need the `gzip` command-line program**,
on `PATH` or supplied with `-gzip /path/to/gzip`. Compressed FASTA uses zlib
directly by default. FASTA in the optional `-resume-mode replay` also uses the
external gzip reader. Compressed separate-FASTA output (`-wsf-gz`) uses zlib
directly and needs no external gzip program. All builds link zlib, including
builds used for BLAST input.

| Sequence input | Plain `.dmp` taxonomy | BLAST taxonomy (`-blasttax`) |
| --- | --- | --- |
| `-blastdb` | `blastdbcmd` | `blastdbcmd` and `sqlite3` |
| Plain or gzip FASTA; optional map plain | No external runtime program | `sqlite3` only |
| FASTA with gzip-compressed map | `gzip` | `gzip` and `sqlite3` |

If either `.dmp` file is compressed, also require `gzip` in that column,
including with `-blastdb`. The three binary/SQLite files used by `-blasttax`
remain extracted, uncompressed files.

Zlib is the only third-party linked library. Standard-library headers such as
`<string>` and `<thread>` come with the compiler. A few POSIX
headers provide process management on Linux/macOS. CladeMers does not need
Python, Boost, OpenMP, an NCBI C++ toolkit build, or a package manager to run.
It uses ordinary 64-bit integers, not a compiler-specific 128-bit extension.

`blastdbcmd` is required only for BLAST sequence input. It understands NCBI's BLAST
database format; CladeMers reads its output and performs exact sequence
subtraction itself. It does **not** run BLAST alignments or use BLAST hit limits.
Using the installed reader avoids maintaining a separate BLAST sequence-database
parser. BLAST taxonomy mode reads the small binary name index directly and uses
the SQLite CLI to read the taxonomy tree.

Build from this directory:

```bash
make clean
make
./clademers --help
```

`make clean` removes `clademers`, `clademers.o` if present, and the
`clademers.dSYM` debug-symbol directory if present. It is safe to repeat and
preserves source files, reference databases, and result directories. Run
`make clean` and `make` as separate commands, especially with parallel Make.

Or use the compiler directly; Make and CMake are not required:

```bash
c++ -O3 -std=c++17 -pthread clademers.cpp -lz -o clademers
```

For a user-local installation, `make install PREFIX="$HOME/.local"` puts the
binary in `$HOME/.local/bin`. Add that directory to your `PATH` if needed.
The same source supports native compilation for **AMD64/x86-64** (Intel or AMD)
and **ARM64/AArch64**. No CPU-specific intrinsics or architecture flags are
required; the selected compiler determines the executable's target. Use
`make CXX=g++` or `make CXX=clang++` to select a compiler. The default flags do
not use `-march=native`, which could make a binary unsuitable for other nodes.
Rebuild after switching the compiler, target architecture, or compiler flags.

Compile on an HPC build/login node or allocated compute node matching the
intended compute architecture and software environment. A macOS executable
cannot be used on Linux; x86-64 and ARM64 also need separate executables.
`blastdbcmd` must match the operating system and architecture too. NCBI supplies
both `x64-linux` and `aarch64-linux` packages in its
[BLAST+ downloads](https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/).
Prefer your HPC's BLAST+ module if available.

Linux/macOS are the intended platforms; native Windows process support and
32-bit ARM validation are outside the current scope. See `VALIDATION.md` for
the platforms actually tested and `HPC.md` for transfer, build, and run commands.

## Run against NT

If you already have BLAST's extracted taxonomy files, use their directory:

```bash
./clademers \
  -blasttax /path/to/blast-taxonomy \
  -taxid 10226 -k 42 -t 30 \
  -blastdb /path/to/nt \
  -out placozoa_k42
```

That directory must contain `taxonomy4blast.sqlite3`, `taxdb.bti`, and `taxdb.btd`
from NCBI's `taxdb.tar.gz`. No `.dmp` download or conversion is needed. The
directory can be the same one that holds NT. `-blasttax` selects these exact
files; it does not change `BLASTDB` or depend on BLAST's taxonomy search path.

Alternatively, the original `.dmp` interface is unchanged:

```bash
./clademers \
  -taxnodes /path/to/taxonomy/nodes.dmp \
  -taxnames /path/to/taxonomy/names.dmp \
  -taxid 10226 \
  -k 42 \
  -t 30 \
  -blastdb /path/to/nt \
  -out placozoa_k42
```

The two dump files can be compressed independently: use
`-taxnodes nodes.dmp.gz -taxnames names.dmp.gz`, or mix plain and gzip files.
This works with BLAST sequence input as well as FASTA. Supply individual dump
files, not the containing `taxdump.tar.gz` archive. Gzip is detected by contents;
no decompressed copy is written to disk. Each dump is streamed once for parsing
and again for content verification before success, with both reader exit
statuses checked. Corruption or a changed dump invalidates the run.

Use **10226 for Placozoa**. The program adds descendants from the selected tree.
Taxid 410657 is the ecological-metagenome group used in your later SRA
metadata filtering; it does not select the placozoan reference lineage.

| Option | Meaning |
| --- | --- |
| `-blasttax` | Directory with the three BLAST taxonomy files; alternative to the two `.dmp` options |
| `-taxnodes` | NCBI `nodes.dmp` or `nodes.dmp.gz`; required in `.dmp` mode |
| `-taxnames` | NCBI `names.dmp` or `names.dmp.gz`; required in `.dmp` mode; supplies scientific names for the target subtree |
| `-taxid` | One target lineage taxid; required |
| `-blastdb` | One complete nucleotide BLAST database prefix; alternative to `-fasta` |
| `-fasta` | Plain or gzip-compressed nucleotide FASTA; alternative to `-blastdb` |
| `-taxidmap` | Optional plain/gzip two-column ID-to-taxid(s) map; only with `-fasta` |
| `-k` | K-mer length, 1–64; default 42 |
| `-t` | Background scanning worker threads; default 1 |
| `-rd`, `--report-descendants` | No value; write a descendant summary now, or later with `-results` |
| `-wsf`, `--write-separate-fastas` | No value; write a FASTA for each taxon with directly assigned markers, plus an index |
| `-wsf-gz`, `--write-separate-fastas-gz` | No value; write those separate FASTAs as `.fasta.gz`; implies `-wsf` |
| `-results` | Completed extraction directory for saved-result mode; requires at least one of `-rd`, `-wsf`, or `-wsf-gz`; also accepts `-out` and `-progress-seconds` |
| `-out` | New output directory; default `clademers_<taxid>_k<k>` |
| `-blastdbcmd` | Optional explicit path to the BLAST+ reader |
| `-sqlite3` | Optional explicit path to the SQLite CLI; only used with `-blasttax` |
| `-gzip` | Optional gzip CLI path for compressed maps/dumps and replay-mode FASTA |
| `-resume-mode` | `direct` (default) or `replay`; choose when starting a new run |

Choose exactly one taxonomy mode; mixing `-blasttax` with either `.dmp` option
is an error. Both modes validate parent links, duplicate taxids, and cycles,
and include internal nodes as well as leaf nodes in the target subtree.

BLAST taxonomy mode reads `TaxidInfo(taxid,parent)` from SQLite and names from
`taxdb.bti`/`taxdb.btd`. BLAST's name database can omit a taxid even when the
tree contains it. In that case, the display label becomes `taxid <number>`
and `target_name_source=taxid_only` is recorded. This does not affect clade
membership or marker sequences. For example, the NCBI snapshot tested here
omitted the Placozoa name at taxid 10226. Missing tree nodes are still errors.

SQLite is opened read-only, its integrity is checked, and the full tree is
validated before sequence extraction. Use an inactive snapshot: nonempty
SQLite `-wal` or `-journal` sidecars are rejected. In both modes, taxonomy files
are fingerprinted on input and rechecked before accepting the result. Dump
fingerprints describe decompressed content; binary/SQLite fingerprints describe
the original file bytes. Plain dump fingerprints retain their previous meaning.

Use taxonomy files matching the NT snapshot when possible. Missing, deleted,
or merged taxids are not automatically repaired; `merged.dmp` is not read.
A numeric reference taxid outside the target subtree, including one absent
from the selected tree, is treated as background. Malformed taxonomy files cause an
error. A target lineage with no descendants is valid.

For an **NT-wide** specificity claim, supply the complete NT prefix, not one
numbered volume, a taxonomic subset, or `core_nt`. CladeMers can also process
other nucleotide BLAST databases; it cannot establish that the input is a
complete or correctly labeled NT release. A single complete alias database
is supported. Lists of multiple database prefixes are not accepted by this
interface; make an alias database first if needed.

Quote paths with spaces normally in your shell. CladeMers supplies BLAST's
additional internal database-name quoting. Literal double quotes in the
database prefix are unsupported. No shell interprets subprocess arguments.

No network connection is used. The output directory must not already exist.

## Run against FASTA

Supply **both target and background references** in one FASTA. Taxids attach
each record to the tree; the tree is still needed to find all descendants.
Choose exactly one of `-fasta` and `-blastdb`.

With taxids in headers:

```bash
./clademers -fasta references.fasta \
  -taxnodes /path/to/taxonomy/nodes.dmp -taxnames /path/to/taxonomy/names.dmp \
  -taxid 10226 -k 42 -t 30 -out placozoa_fasta_k42
```

Accepted header examples (keys are case-sensitive; no spaces inside a taxid list):

```text
>ref1 taxid=10227
>ref2 [taxid=10227]
>ref3 taxids=10227,9606
>ref4 [taxids=10227;9000001]
>ref5|kraken:taxid|10227
>kraken:taxid|10227|ref6
```

`taxid=` and `taxids=` must be separate whitespace-delimited fields after the
first identifier; brackets are optional. Kraken annotations belong in the
identifier and contain one taxid. If a header has multiple declarations, their
normalized taxid sets must agree. Malformed recognized annotations are errors.
CladeMers does not infer taxids from bare numbers, organism names, accessions,
or arbitrary header formats.

For a separate mapping file, add `-taxidmap sequence_taxids.tsv` or
`-taxidmap sequence_taxids.tsv.gz`:

```bash
./clademers -fasta references.fasta -taxidmap sequence_taxids.tsv.gz \
  -blasttax /path/to/blast-taxonomy \
  -taxid 10226 -k 42 -t 30 -out placozoa_mapped_k42
```

The FASTA and mapping file can be compressed independently, and either example
can use either taxonomy option. The decompressed map has exactly two
whitespace-delimited columns, with an optional `sequence_id taxid` or
`sequence_id taxids` first row:

```text
sequence_id     taxids
ref1            10227
ref2            10227;9000001
ref3            10227
ref3            9606
ref4            0
```

- The key is the **exact first whitespace-delimited FASTA identifier**.
  Accession versions and pipe-delimited identifiers are not stripped.
- Multiple rows for the same ID are combined as a union. Above, `ref3` is
  mixed target/background and is subtracted. Lists accept commas or semicolons.
- Blank lines and full lines starting with `#` are ignored. Inline comments
  and extra columns are not supported. Unused map IDs are allowed.
- The selected map is authoritative: an ID missing from it is unknown/background,
  even if its header has a taxid. When a mapped ID also has an explicit header
  annotation, the sets must agree or the run fails.
- Missing header assignments, `0`, and `N/A` are unknown/background. Any record
  containing a non-target or unknown assignment is background, including a
  mixed target/background record. Positive taxids absent from the supplied
  tree are also background. Nothing is silently discarded from the scan.
- Standard NCBI four-column accession-to-taxid files are **not** this map format.
  Convert them to two columns matching your FASTA's exact IDs first. An ordinary
  NCBI FASTA header often needs this external mapping; names are not sufficient.

The FASTA must be a **regular file**, plain or gzip-compressed, not stdin or a
pipe. Pass `.fasta.gz` or `.fa.gz` directly with the same `-fasta` option:

```bash
./clademers -fasta references.fasta.gz -taxidmap sequence_taxids.tsv.gz \
  -taxnodes nodes.dmp.gz -taxnames names.dmp.gz \
  -taxid 10226 -k 42 -t 30 -out placozoa_gzip_k42
```

Omit `-taxidmap` to use header taxids. All four text inputs must be regular
files and can be compressed independently. Compression is detected from each file's gzip
header bytes, so other filename extensions also work. A `.gz` file without a
gzip header is rejected. Concatenated gzip members are supported. The readers
stream directly into the parser: **no uncompressed copy is written to
disk**, and the compressed file is preserved. CRC errors, truncated streams,
and nonzero reader exit statuses invalidate the run. A compressed map uses the
same two-column text format, comments, header, duplicate-ID union rules, and
exact ID matching as a plain map.

Wrapped sequences, blank lines, CRLF, lowercase, and a missing
final newline are supported. Spaces/tabs in sequences are removed; IUPAC
ambiguities, `U`, and `-` interrupt ACGT windows and count toward coordinates.
Other sequence characters, empty records, and control characters in headers
(including merged defline separators) are rejected. Convert merged taxonomy
annotations into explicit taxid lists; do not discard their non-target labels.

The file is scanned twice, decompressing it on each pass if needed, and the
exact **decompressed content bytes** are fingerprinted before stripping whitespace.
Different fingerprints invalidate the result. For plain input this is the
original file's byte fingerprint, as in 0.3.0. Gzip header timestamps and
compression settings are not part of the content fingerprint.
The map is loaded in memory, then streamed again before completion to verify
its decompressed content fingerprint and compression mode. The verification
does not build a second map. A compressed map is therefore decompressed twice;
both reader invocations must succeed. Plain-map fingerprints are unchanged.
Keep all input files unchanged during a run.
The source record number is retained, so duplicate FASTA IDs are supported.
Specificity covers only the references you supply; a target-only FASTA cannot
establish absence from other lineages.

## Exactly what is computed

1. Read the taxonomy tree and identify the target taxid and all descendants.
2. Read target-associated references. For BLAST input, ask `blastdbcmd` for
   that explicit list using `-no_taxid_expansion`. For FASTA, scan the full
   file and classify each record's assigned taxids against the same tree.
3. Enumerate every A/C/G/T k-mer in those target references. Store the complete
   sequence identity and one representative source occurrence in memory.
   Always accumulate the LCA across all of that k-mer's taxonomic assignments;
   later occurrences can broaden its label.
4. Stream **all** input records and remove any candidate that occurs in
   background, including unknown or mixed-taxonomy records.
5. Check BLAST producer exit status and database metadata, or FASTA pass
   fingerprints and the mapping file; recheck taxonomy files, then write
   the retained set.

The set operation is:

```text
retained = distinct canonical k-mers in target references
           minus all k-mers seen in background references
```

**Canonical** means the lexicographically smaller of the forward sequence and
its reverse complement. Either strand in background excludes the candidate.
Lowercase A/C/G/T is accepted. Other bases interrupt the window; windows never
cross record boundaries. Parallel chunks overlap by `k-1` bases so a window
crossing a chunk boundary is still checked.

BLAST export uses `%a`, `%T`, `%X`, and `%s`: accession, principal taxid, leaf
taxids, and unwrapped sequence. **All deflines are exported**, without
`-target_only`, so a merged sequence with a non-target defline is subtracted.
The full background export also includes `%o` (OID) for progress accounting;
it does not change which rows are scanned.
A row is exempt from subtraction only when its principal taxonomy and every
reported leaf taxid lie inside the target clade. An absent optional leaf list
(`N/A`, `0`, or empty) is allowed when the principal annotation is in the clade.
An unknown principal taxid or a mixed leaf list is background. This policy is
conservative: an ancestral principal annotation with target-only leaves is
also background.

All `2*k` bits determine identity: **42 bases use 84 bits**, not a 64-bit hash.
A small immutable prefilter speeds up negative checks; only full sequence
equality can exclude a marker. Counts, abundance cutoffs, low-complexity
filters, and sampling do not alter the returned set. Singleton markers and
low-complexity sequences are retained.

## Progress and estimated time remaining

Progress is enabled automatically on **stderr**, so it appears in the existing
Slurm log. Updates include percent complete, elapsed time, approximate ETA, and
processed/total counts for the **named operation**. For example (illustrative):

```text
      [progress] BLAST background scan: 40.0% | elapsed 00:20:00 | ETA ~00:30:00 | 400000000/1000000000 sequence bases
```

The default interval is 30 seconds. Add `-progress-seconds 10` for more frequent
updates (integer seconds, 1–3600). Start/completion lines are always printed.
Updates happen at processing boundaries, so a single long record, sorting,
or a stalled reader can delay the next line. For an interactive run, append
`2> clademers.log` to the command and use `tail -f clademers.log` to watch it.

- **BLAST background scan:** uses total bases reported by `blastdbcmd -info`.
  OIDs prevent consecutive deflines of the same sequence from inflating progress;
  all deflines still participate in subtraction. Reported BLAST nucleotide
  lengths can be approximate. Missing/unrecognized totals, counts exceeding the
  total, or unexpected OID ordering disable the estimate without stopping the scan.
- **Plain text inputs:** use file bytes. The second FASTA pass and map/dump
  rechecks use content sizes measured during the first pass, including gzip inputs.
- **First gzip passes, BLAST target extraction, SQLite export, and sorting:**
  a reliable total is not available in advance, so the log explicitly says
  `percent unavailable` and `ETA unavailable`. No extra counting/decompression
  pass or temporary uncompressed file is introduced. Gzip footer sizes are not
  used; concatenated streams and files over 4 GiB do not provide a safe total there.
- **Collecting/writing markers:** use candidate/retained-marker counts.

ETA extrapolates the operation's average rate and can change with sequence lengths,
CPU load, and storage throughput. It does **not** estimate the remaining whole run.
The FASTA background-input percentage covers reading/queueing; workers and input
snapshot verification finish afterward. The BLAST scan reaches 100% only after
its reader succeeds and workers finish. A completed read/recheck line describes
I/O completion; fingerprint comparisons can still reject changed inputs.
Only `SUCCESS` plus a complete manifest establish successful completion of the run.

## Descendant-specific markers

Every new run saves LCA labels automatically. Add `-rd` (or
`--report-descendants`, with no value) to also write the summary immediately:

```bash
./clademers -blastdb /path/to/nt -blasttax /path/to/blast-taxonomy \
  -taxid 10226 -k 42 -t 30 -rd -out placozoa_descendants_k42
```

This works with FASTA, header or mapping-file taxids, both taxonomy sources,
and compressed text inputs. Marker sequences, IDs, and representative source
coordinates are unchanged. **Every run**, even without `-rd`, writes these
columns at the end of `markers.tsv`:

| Column | Meaning |
| --- | --- |
| `specificity_taxid` | Lowest common ancestor (LCA) of **all observed taxonomic assignments** for that retained k-mer |
| `specificity_name` | Scientific name of that clade, or `taxid <number>` when a descendant name is missing |

FASTA descriptions also save the LCA taxid, using the existing `taxid=` header
convention. For example:

```text
>km_10226_k42_1 taxid=10227
GATCCGAATGCGTACCTGATCGATGCTAGTCCGATAGCTTAC
```

Here `10226` in the marker ID is the requested target; `taxid=10227` is this
marker's most-specific assignment. The first whitespace-delimited identifier
stays unchanged. The saved taxid is the LCA, not a list of all source taxids.

Each k-mer is reported once, under its most specific supported clade. A marker
observed only in species A is assigned to A. If it also occurs in sibling
species B, it is assigned to their common ancestor. A marker shared across the
selected target's main branches may be assigned to the target itself. The
assignment includes every occurrence, reverse complements, multiple taxids on
a record, and BLAST principal **and** optional leaf taxids. A broad principal
annotation limits the assignment even when a narrower leaf annotation is also
present. The representative source accession alone does not determine the label.

`descendant_summary.tsv` is sorted by taxid and contains every taxon in the
selected subtree, including the requested target:

| Column | Meaning |
| --- | --- |
| `taxid`, `parent_taxid`, `scientific_name` | Taxonomy node and its parent in the supplied tree |
| `assigned_kmers` | Number of markers whose most-specific assignment is exactly this node |
| `assigned_kmers_incl_descendants` | Number assigned to this node **or any of its descendants** |

The sum of `assigned_kmers` equals the number of retained markers. Inclusive
`assigned_kmers_incl_descendants` counts overlap between ancestors and descendants; the
target's inclusive count equals the total. Zero counts do not establish that a
taxon is absent from the reference. To extract a node's compact marker list,
filter `markers.tsv` on `specificity_taxid`.

All markers still pass the full background exclusion check. Unknown, mixed,
or out-of-target assignments remove matching candidates. Specificity is relative
to the supplied reference snapshot and its annotations: it does **not** mean a
marker occurs in every species or genome inside its assigned clade. Broad or
missing annotations can limit resolution or exclude otherwise useful markers.

Assignments accumulate during the existing target pass, with no additional NT
scan. Every run retains the selected subtree's parent/depth/name data and one
taxid per candidate. Checkpoints always preserve these assignments.

### Generate the report later

After a run completes, generate the same summary using only its saved results:

```bash
./clademers -rd -results placozoa_k42 -out placozoa_descendant_report
```

Omit `-out` to use `placozoa_k42/descendant_report`. The report directory must
not already exist. This operation leaves the extraction files unchanged.
With `-rd`, it writes `descendant_summary.tsv`, a report `manifest.tsv`, and its
own `SUCCESS`. With `-wsf` or `-wsf-gz`, it also (or instead) writes the separate FASTAs and
index; metadata and `SUCCESS` are written in either mode.
It produces the same summary bytes as `-rd` on the original extraction.

Keep these four files together (the directory can be moved to another computer):

- `markers.fasta`: canonical marker sequences and saved LCA taxids.
- `target_taxonomy.tsv`: the selected subtree's taxids, parents, and scientific names.
- `manifest.tsv`: completion status, counts, format version, and file fingerprints.
- `SUCCESS`: the completed-extraction marker.

Reporting does **not** read NT, source FASTA, mapping files, taxonomy dumps,
BLAST taxonomy, `markers.tsv`, or checkpoint files. It needs no BLAST+, SQLite,
or gzip executable. Keep the saved marker FASTA uncompressed and unedited;
reporting validates its fingerprint and record count, plus the taxonomy snapshot.
The fingerprints detect accidental changes, not malicious tampering.

Labels and counts refer to the original reference/taxonomy snapshot. A saved
LCA cannot reconstruct all original occurrence taxids under a different tree.
This portable results format starts with **0.11.0**. Earlier outputs lack the
required catalog metadata/snapshot and cannot use `-results`; this release
requires a new extraction for them. Report-only work is not checkpointed; if
interrupted, rerun it into a new report directory. NT is never rescanned by this
mode.

See [the synthetic descendant demo](examples/descendant_demo/README.md) for a
small runnable example and its expected report.

### Write separate FASTAs per taxon

Use `-wsf` (or `--write-separate-fastas`, with no value) during extraction, or
export them later from any completed **0.11.0 or newer** extraction:

```bash
./clademers -wsf -results placozoa_k42 -out placozoa_separate_fastas
# Add -rd if you also want the descendant summary:
./clademers -wsf -rd -results placozoa_k42 -out placozoa_report_with_fastas
# Or compress the separate FASTAs (no need to repeat -wsf):
./clademers -wsf-gz -rd -results placozoa_k42 -out placozoa_report_with_gzip_fastas
```

Each marker goes into exactly one file, matching its saved most-specific LCA.
For taxid 10227, that file is `taxon_fastas/taxid_10227.fasta`. Headers, marker
IDs, sequences and within-file sequence order match the combined catalog.
Markers assigned to child taxa are not copied into ancestor files. A marker
assigned to the requested target itself goes into that target's file.

`-wsf-gz` / `--write-separate-fastas-gz` writes
`taxon_fastas/taxid_10227.fasta.gz` instead. It implies `-wsf`; if both flags
are present, compression is enabled regardless of their order. Use it during
extraction, with `-results`, or when resuming an incomplete run. Compression
uses zlib at level 6, with no new dependency or external gzip process. The
decompressed FASTA content is identical to the plain separate output. The
combined `markers.fasta`, `markers.tsv`, and report/index files remain plain.

The index `taxon_fastas.tsv` contains `taxid`, `scientific_name`,
`assigned_kmers`, and `fasta_file` (relative to the output directory), sorted by
taxid. Compressed exports use `.fasta.gz` paths in this index. Each row
corresponds to a nonempty FASTA. Its count matches
`assigned_kmers` in the descendant summary, not
`assigned_kmers_incl_descendants`. Taxa with zero direct assignments have no
FASTA/index row; an empty catalog produces an empty `taxon_fastas/` directory
and a header-only index. Numeric filenames avoid scientific-name collisions
and punctuation/path issues.

`-wsf` and `-rd` are independent. The combined `markers.fasta` and `markers.tsv`
remain part of every extraction. The separate files add another copy of the
marker FASTA content plus filesystem overhead, and a large clade may create
many files. The writer keeps at most 32 taxon files open at once. If a compressed
file is closed and reopened, another standard gzip member is appended. Use a
reader that handles concatenated gzip members, such as `gzip -dc`, Python's
`gzip`, or zlib's gzip reader. Many small members can reduce compression
efficiency; tiny FASTAs can become larger when compressed. Saved-result export
streams the marker catalog once and never reads NT or original inputs.

Separate FASTAs are staged in `taxon_fastas.partial/` until all records and
fingerprints have passed validation. As with other outputs, use only results
with `SUCCESS` and a complete manifest. Output-phase checkpoint recovery
regenerates the files; stale generated files are replaced without appending
duplicate markers. Unexpected contents in a generated directory are rejected.
Both the separate-FASTA choice and its compression setting are saved in new
checkpoints. Adding `-wsf-gz` with `-resume` switches an incomplete plain run to
compressed output; repeating `-wsf` does not disable saved compression.
Output recovery removes stale generated files in either format before publishing
the new export. Compressed writes and finalization are checked for errors.

## Checkpoint and resume

Checkpoints are enabled by default, with a **300-second interval** plus saves at
phase boundaries. Keep the output directory on persistent storage accessible
from subsequent jobs. After interruption, use the **same CladeMers release**:

```bash
./clademers -resume /scratch/yourname/placozoa_k42 -t 30
```

The checkpoint restores the sequence/taxonomy paths, taxid, k, candidates,
representative source coordinates, output options and LCA assignments,
processed-record counts, and exclusion flags.
Do not repeat `-out`, input arguments, `-taxid`, or `-k` with `-resume`. You may
add `-rd` for the summary and `-wsf` or `-wsf-gz` for separate FASTAs on an
incomplete run; LCA data was already saved.
You may also override `-t`, `-progress-seconds`, `-checkpoint-seconds`, and reader executable
paths (`-blastdbcmd`, `-sqlite3`, `-gzip`). Taxonomy and mappings are reloaded and
compared to the saved fingerprints. BLAST metadata and reader version must match.

**Direct resume is the default**, in both target loading and
background scanning:

| Sequence input | How a restart reaches saved work |
| --- | --- |
| NT / BLAST | A temporary OID-range alias asks `blastdbcmd` to begin at the last saved sequence OID. Only the saved defline rows at that OID are re-exported and verified; earlier sequences are not exported. Original alias filters and volume ordering are retained. |
| Plain FASTA | Seeks to the saved byte offset and restores the next header already read by the parser. |
| Gzip FASTA | Seeks to a saved DEFLATE block boundary, restores up to 32 KiB of dictionary plus pending bits and CRC state, then decompresses locally to the saved record boundary. |

Gzip access points are built during the normal scan and stored inside the
checkpoint. There is no separate indexing pass, persistent index to manage, or
uncompressed FASTA copy. Points are sampled at block boundaries about 1 MiB of
content apart; the replay distance also depends on the compressor's block sizes.
An unusually long DEFLATE block can require a longer local replay. Concatenated
members work, and resumed reading still checks member CRCs and lengths. Log
messages report the selected BLAST OID, plain byte offset, or compressed byte
offset and local replay distance.

Direct resume requires **unchanged reference files at their original paths**.
It compares canonical paths, sizes, and nanosecond modification/change timestamps
for the FASTA or all resolved BLAST alias, volume, and filter files. Do not update,
replace, touch, move, or recopy these files between jobs. These metadata checks
detect ordinary changes without rereading NT; they are not full content checksums
or protection against deliberate timestamp manipulation. Keep reference releases
immutable and retain their checksum manifests. BLAST prefixes are resolved to
absolute paths at startup; use an absolute complete NT prefix for clarity.

For an unusual database alias unsupported by direct file discovery, start a
**new** run with `-resume-mode replay`. This retains the older prefix-replay
behavior and exact saved-prefix verification, including external gzip for FASTA.
The mode is saved in the checkpoint and cannot be changed on resume. In replay
mode, bare BLAST names also require the same working directory and `BLASTDB` value.

Taxonomy/map loading and checkpoint-table loading still take time on restart.
Completed target collection is skipped. A checkpoint made after input verification
restarts sorting/output with input checks repeated; partial marker files are
rewritten. Direct mode does not reread the FASTA sequence file in this output phase.
**Older checkpoints require their original release to resume.** Version 0.13.0
uses checkpoint state format 5 to save the separate-FASTA compression setting.
Completed 0.11.0 and newer extractions work with `-results` and any combination
of `-rd`, `-wsf`, and `-wsf-gz`; they do not need new extraction or checkpoint
conversion.

`-checkpoint-seconds 600` changes the interval; `0` disables new saves.
An interval is a minimum between saves, not a guaranteed recovery-time bound:
saves occur at complete-record boundaries, after queued workers drain. Large
records, checkpoint writes, and sorting can delay signal handling. Taxonomy/map
setup is repeated on restart; an interruption before the initial checkpoint
requires restarting with a new output directory. Runs begun with older releases
without checkpoints cannot be resumed retroactively.

Send **SIGUSR1**, **SIGTERM**, or **SIGINT** to the CladeMers process to request a
safe save and pause (exit status **75**, `manifest.tsv` status `paused`). SIGKILL,
node loss, or a hard scheduler kill cannot trigger a final save: resume uses the
last committed checkpoint. No automatic job submission/requeue is performed.
See [the Slurm checkpoint example](examples/nt_checkpoint.slurm) for an advance
signal, forwarding it to CladeMers, waiting for the save, and resuming on the next
submission. Set its warning lead time longer than your observed record/checkpoint
latency and keep the result path constant across submissions.

Checkpoint files use explicit byte order, so they do not depend on AMD64 versus
ARM64 native struct layout. Compile CladeMers and install its readers for the
new node's architecture; do not copy an incompatible executable. The input paths,
reader output/version, and release must still match.

Writes use a temporary file, file sync, atomic rename, and directory sync. A
failed/unfinished save leaves the previous committed state usable; abandoned
`.partial` files are ignored. Checksums detect accidental checkpoint damage.
An advisory lock prevents concurrent writers. Use a filesystem supporting those
operations, and retain **both** checkpoint files and the result directory.
During target loading, each save writes the growing table. Afterward, only the
compact exclusion flags and metadata change. Budget roughly **56 bytes per
candidate plus source strings**, with room for a second table during atomic
replacement/phase transition, in addition to result files. This disk estimate
does not describe RAM requirements; abandoned partial files can use additional
space after an interruption. Checkpoints remain after success;
`make clean` preserves them. A completed directory cannot be resumed.

## Outputs

Treat a run as complete only when **`SUCCESS` exists and `manifest.tsv` says
`status = complete`**. An interrupted or failed run is not a specificity result.

| File | Contents |
| --- | --- |
| `markers.fasta` | All retained canonical k-mers; each header includes `taxid=<LCA>` |
| `markers.tsv` | Marker sequence, one source accession/ID, taxids, coordinate, strand; FASTA also includes record number; always includes specificity taxid/name |
| `descendant_summary.tsv` | With `-rd`: counts for every taxon in the selected subtree, including the target and zero-count taxa |
| `taxon_fastas/taxid_<id>.fasta` or `.fasta.gz` | With `-wsf` or `-wsf-gz`: markers assigned directly to this LCA taxid; nonempty taxa only |
| `taxon_fastas.tsv` | With either separate-FASTA option: taxid, scientific name, direct count, and relative FASTA filename |
| `target_taxids.txt` | Sorted target subtree used for extraction |
| `target_taxonomy.tsv` | Sorted target subtree with parent taxids and scientific names, saved for later reporting |
| `manifest.tsv` | Options, versions, taxonomy fingerprints, counts, timing, completion status |
| `database_info.txt` | BLAST sequence input only: database description, timestamp, volumes, and sizes |
| `commands.log` | Exact subprocess arguments; absent when no subprocess is used |
| `blastdbcmd.log` | BLAST sequence input only: BLAST diagnostics |
| `taxonomy.log` | SQLite diagnostics; present in BLAST taxonomy mode |
| `gzip.log` | External gzip diagnostics, when used for maps/dumps or replay-mode FASTA |
| `resume_range.nal` | Generated BLAST OID-range alias during direct resume |
| `checkpoint/state.bin` | Latest committed restart state; automatically saved unless disabled |
| `checkpoint/candidates.bin` | Immutable candidate/provenance table after target loading completes |
| `run.lock` | Advisory writer lock; automatically released on process exit, including forced termination |
| `SUCCESS` | Written atomically only after successful processing and output |

The manifest records `sequence_source` (`blast` or `fasta`),
`taxonomy_source` (`dmp` or `blast`), input file paths,
byte sizes, fingerprints, the target label's source, and the SQLite CLI version
when used. It also records `progress_seconds`, `progress_scope`, and
`blast_progress_total_bases` when a BLAST size hint is available.
Portable result metadata includes `results_format=1`, `lca_saved=yes`, and
`markers_fasta_*` / `target_taxonomy_*` byte counts and FNV-1a64 fingerprints.
Restart metadata includes `checkpoint_seconds`, `checkpoint_phase`,
`checkpoint_records`, `resume_count`, `resume_mode`, and `session_started_utc`.
Direct runs also describe `sequence_snapshot_check`. `elapsed_seconds`
includes time recorded in previous checkpoints plus the current session; work
lost after a hard kill is not included in that saved-time total.
FASTA runs also record the annotation source, input bytes and FNV-1a64
fingerprint, matching passes, record count, and optional map fingerprint/ID count.
`fasta_compression` is `none` or `gzip`; `fasta_fingerprint_scope` is
`decompressed_content`. Both `fasta_bytes` and `fasta_fnv1a64` describe the
decompressed content, not the gzip container. `fasta_reader` distinguishes
`plain_file`, `indexed_zlib`, and `gzip_command`; `zlib_version` records the linked
library. When external gzip is used, its path and invocations are recorded in
`commands.log` (two per compressed input in an uninterrupted full run).
Map provenance similarly records `taxidmap_compression`,
`taxidmap_fingerprint_scope=decompressed_content`, `taxidmap_bytes`,
`taxidmap_fnv1a64`, and `taxidmap_passes_match=yes` after successful verification.
Map byte counts/fingerprints describe decompressed text, not its gzip container.
Dump provenance similarly records `nodes_compression`, `names_compression`,
`nodes_fingerprint_scope`, and `names_fingerprint_scope`; the scopes are
`decompressed_content`, and their byte counts/fingerprints describe that content.
BLAST taxonomy file fingerprint scopes are `file_bytes`.
`taxonomy_inputs_match=yes` records successful final verification in either mode.
`fasta_unassigned_taxid_records` counts records with missing assignments or
`0`/`N/A`; it does not count other positive taxids outside the target subtree.
The existing `*_defline_rows` count keys count FASTA records in FASTA mode.
Fingerprints are not cryptographic integrity guarantees.

The TSV's `start_1based` is the leftmost position of the reference window.
`strand = -` means reverse-complement that window to obtain the emitted marker.
The source is the first target occurrence encountered, **not all occurrences**.
For FASTA, `source_accession` is the complete first header token,
`source_taxids` is the normalized assigned list, `source_leaf_taxids` is `N/A`,
and a trailing `source_record_1based` column identifies the original record.
The seven-column BLAST TSV format is unchanged. Coordinates refer to the
sequence after stripping whitespace, including ambiguity/gap characters.
FASTA and TSV are sorted by sequence and deterministic across worker counts.
Sequential marker IDs apply to that particular catalog, not across NT updates.

If every candidate is shared, empty FASTA and header-only TSV files are a
valid result. If there are no target k-mers, or no background records at least
`k` bases long, the program stops with an error rather than an uninformative
specificity claim. Failed jobs keep diagnostics and any partial files.

## Memory, threads, and HPC use

- **Every distinct target k-mer stays in RAM** until the run ends, including
  candidates excluded during the background scan. Broad target lineages can
  be too large. There is no disk-backed candidate table or memory limit yet.
- With BLAST input, the full NT export is streamed through a pipe, not written to disk. Only
  the current exported row is buffered by the producer, plus a bounded queue
  of approximately 1 MiB chunks and worker buffers. A very long reference
  sequence still needs a correspondingly large input-line buffer.
- FASTA input buffers one complete record plus its current line. It reads the
  entire file twice; gzip input is decompressed on each pass without disk staging.
  A supplied taxid map is held fully in RAM, including unused IDs, and streamed
  again for verification; compressed maps are decompressed on both reads.
  Compression saves disk space, not the RAM used by loaded map entries.
  Large accession maps can therefore require substantial extra memory.
- Taxonomy dumps are also read twice: once to load/validate the tree and target
  name, and once to recheck content. The second read uses a line buffer without
  rebuilding the tree. Gzip saves disk space, not taxonomy-tree RAM.
- `-t 30` starts 30 scanning workers. The producer thread and `blastdbcmd`
  process can also consume CPU; an allocation of 32 CPUs accommodates these.
  Indexed FASTA decompression runs in the producer thread. Compressed maps/dumps
  and replay-mode FASTA use a gzip reader process.
  Decompression, target loading, input parsing,
  database export, sorting, and output remain serial. The
  database reader or storage bandwidth may limit scaling.
- Candidate memory is not multiplied by the thread count. Workers share a
  read-only candidate table and set atomic exclusion flags.
- Output can be large: both full FASTA and TSV are written without compression.
  Compress completed outputs separately if desired.
- Direct checkpoint resume skips the completed sequence prefix, with local gzip
  block replay or BLAST OID-group verification. Setup and checkpoint loading
  still take time. Job resubmission is manual.

Keep NT unchanged throughout the run. Direct mode checks reference-file metadata
as well as `blastdbcmd -info`; neither is a full database checksum. Record the
NT release and retain the supplied database's checksum manifest separately.
The FNV-1a64 fingerprints in `manifest.tsv` identify binary/SQLite taxonomy-file
bytes and decompressed FASTA/map/dump content for
accidental-change checking; they are not cryptographic hashes.

See `examples/nt.slurm` for `.dmp` inputs or
`examples/nt_blast_taxonomy.slurm` for BLAST taxonomy inputs. Their RAM and
time requests are placeholders, not measurements of a full NT run.

## Try a tiny example first

The example sequences are synthetic test data, not placozoan markers. The
following two commands need neither BLAST+ nor SQLite and each retain one of
two candidate 42-mers (use fresh output directories):

```bash
./clademers -fasta examples/references_taxids.fasta \
  -taxnodes examples/nodes.dmp -taxnames examples/names.dmp \
  -taxid 10226 -k 42 -out demo_headers
./clademers -fasta examples/references.fasta -taxidmap examples/taxid_map.tsv \
  -taxnodes examples/nodes.dmp -taxnames examples/names.dmp \
  -taxid 10226 -k 42 -out demo_mapping
```

To compare with BLAST input, with BLAST+ on `PATH`:

```bash
demo_dir=$(mktemp -d /tmp/clademers-demo.XXXXXX)
makeblastdb -in examples/references.fasta -dbtype nucl \
  -parse_seqids -taxid_map examples/taxid_map.tsv -blastdb_version 5 \
  -out "$demo_dir/demo_db"
./clademers -taxnodes examples/nodes.dmp -taxnames examples/names.dmp \
  -taxid 10226 -k 42 -t 2 -blastdb "$demo_dir/demo_db" -out "$demo_dir/result"
cat "$demo_dir/result/markers.tsv"
```

## Developer tests

Tests use Python 3's standard library; Python is not needed for normal runs.

```bash
make test
make test-blast BLAST_BIN=/path/to/ncbi-blast/bin
# Optional: also exercise an existing, genuine NCBI BLAST taxonomy snapshot:
make test-blast BLAST_BIN=/path/to/ncbi-blast/bin BLAST_TAXDB=/path/to/blast-taxonomy
```

The second command includes an actual version-5 database built by `makeblastdb`.
The suite compares results against a direct DNA-string oracle, exercises bit
boundaries (including 32, 33, 42, and 64 bases), checks thread determinism,
coordinates, reverse complements, chunk boundaries, unknown/mixed taxonomy,
empty results, producer failures, malformed files, and stale output protection.
FASTA tests cover headers, map unions/conflicts, exact ID matching, missing
assignments, wrapping, CRLF, malformed input, duplicate-ID provenance, and
operation without BLAST+. With `gzip` available, tests compare all 16 plain/gzip
FASTA/map/nodes/names combinations, compressed dumps with genuine BLAST input,
concatenated members, large lines and buffer boundaries,
missing final newlines, corrupt/truncated data, reader failures, and changed
second-pass content, including map/dump changes detected during final verification.
`CLADEMERS_GZIP` selects a particular gzip CLI for developer tests.
With `sqlite3` available, it also compares the two taxonomy modes, validates
the binary name format, and checks corrupt/missing files, invalid trees,
missing internal names, failed SQLite readers, and changed taxonomy inputs.
Tests needing unavailable tools/data are explicitly skipped. `CLADEMERS_SQLITE3`
can select a particular SQLite CLI for developer tests. No test downloads NT.
Restart tests send real SIGUSR1/SIGTERM/SIGINT/SIGKILL signals, compare resumed
markers/provenance/counts with uninterrupted runs, and reject changed inputs,
corrupt state, and concurrent writers. `CLADEMERS_RESUME_BINARY` optionally
selects a second build to test checkpoint exchange across architectures.
Direct-resume tests check plain byte seeking, indexed gzip restart in both phases,
gzip optional headers/CRC/truncation, and real BLAST OID ranges across nested
filtered aliases with repeated deflines. Older reader compatibility tests select
`-resume-mode replay` explicitly.
Descendant tests compare all-occurrence LCAs with an independent string/set
oracle and separate per-descendant runs. They cover sibling sharing, broad and
mixed annotations, both taxonomy sources, compressed inputs, missing names,
empty reports, thread determinism, and checkpoint recovery of labels/counts.
Deferred-report tests remove original inputs and reader executables, move only
the four required result files, and compare immediate/later summaries. They
reject incomplete/legacy results, damaged catalogs, invalid saved trees, and
unsupported CLI combinations. Separate-FASTA tests check exact partitions,
identical decoded immediate/later exports, numeric filenames, empty output,
checkpoint recovery, and more than 100 taxa under a 64-file-descriptor limit.
Gzip tests validate all concatenated members and checksums, compression restored
from checkpoints, switching formats during output recovery, prior-version
portable exports, and finalization failure without publishing partial results.

For a guided introduction, see [Understanding CladeMers](CODE_WALKTHROUGH.md):
the biological model, a worked example, C++ vocabulary, function-by-function
reference, threading, input handling, and hands-on exercises. That guide describes
the 0.7.0 implementation; use the sections above for current reporting and restart behavior.

To study the implementation, read `clademers.cpp` from top to bottom. It is
organized into annotated sections. `each_kmer()` implements rolling
sequence encoding; `classify_record()` defines the shared taxonomy policy;
`read_fasta()` and `parse_record()` supply records; `run()` puts
the biological workflow together. Most remaining code handles safe process
execution, bounded parallel work, and auditable output.

## Interpretation and scope

"Lineage-specific" here means **present in at least one target reference and
absent from screened background references in this database snapshot**.
It does not mean conserved across the target lineage, absent from unsequenced
organisms, or immune to reference contamination and annotation errors.

This program discovers reference markers. It does not submit PebbleScout
queries, filter SRA metadata, or verify exact matches in reads. Those are
separate downstream steps. A PebbleScout candidate hit still needs exact-read
verification and biological evaluation before being called a placozoan detection.

## BLAST documentation behind the interface

- [Taxonomy filtering and `-no_taxid_expansion`](https://www.ncbi.nlm.nih.gov/books/NBK569839/)
- [Retrieval by taxid and non-redundant deflines](https://www.ncbi.nlm.nih.gov/books/NBK569840/)
- [BLAST database paths containing spaces](https://www.ncbi.nlm.nih.gov/books/NBK569860/)
- [NCBI sequence formatter implementation](https://github.com/ncbi/ncbi-cxx-toolkit-public/blob/main/src/objtools/blast/blastdb_format/seq_formatter.cpp)
- [NCBI BLAST SQLite tree schema and query](https://github.com/ncbi/ncbi-cxx-toolkit-public/blob/main/src/objtools/blast/seqdb_reader/tax4blastsqlite.cpp)
- [NCBI binary taxonomy name format](https://github.com/ncbi/ncbi-cxx-toolkit-public/blob/main/src/objtools/blast/seqdb_reader/seqdbtax.cpp)
- [SQLite command-line reader](https://sqlite.org/cli.html)
- [NCBI alias range implementation](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/doxyhtml/seqdbalias_8cpp_source.html)
- [Zlib manual: block boundaries, dictionaries, and pending bits](https://zlib.net/manual.html)
- [Zlib random-access example and access-point design](https://github.com/madler/zlib/blob/master/examples/zran.c)
