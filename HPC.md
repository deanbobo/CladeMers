# CladeMers on an HPC

The source archive contains C++ source, build instructions, developer tests,
and synthetic examples. It contains no compiled executable, BLAST+ software,
reference database, or actual Placozoa marker catalog.

## Select the build environment

Use a node and compiler toolchain matching your intended compute nodes:

| `uname -m` | Architecture | Native build |
| --- | --- | --- |
| `x86_64` | AMD64/x86-64, used by AMD and Intel CPUs | `make CXX=g++` |
| `aarch64` or `arm64` | 64-bit ARM (ARM64/AArch64) | `make CXX=g++` |

The source is the same for both. The executable is specific to the operating
system, architecture, and runtime libraries it was built for. Use your HPC's
compiler module or ask its documentation which toolchain to load; module names
vary by site. GCC 9+ or a recent Clang with C++17 support is needed. Make builds
the selected compiler's target; it does not select a remote compute-node
architecture automatically. If login and compute architectures differ, build
on an allocated node of the intended architecture or use a suitable cross-toolchain.

Version 0.13.0 also requires **zlib development headers and the library**, matching
the target architecture. Load your site's zlib module/package before building.
The OS development packages are commonly `zlib1g-dev` (Debian/Ubuntu) or
`zlib-devel` (RHEL/Rocky/AlmaLinux). The Makefile adds `-lz`; the shared zlib
runtime must also be available on the compute nodes.

BLAST sequence input requires BLAST+ 2.15 or newer, with a matching `blastdbcmd`
on `PATH`. FASTA sequence input does not need BLAST+.
Use your site's BLAST+ module when available. NCBI's
[BLAST+ download directory](https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/)
provides Linux packages labeled `x64-linux` and `aarch64-linux`; choose the one
matching your compute nodes. Source builds do not require Python. Developer
tests optionally use Python 3's standard library and `makeblastdb`.
If using BLAST taxonomy files with `-blasttax`, also load or install the
`sqlite3` command-line reader. Check it with `sqlite3 -version`; if it is
outside `PATH`, add `-sqlite3 /path/to/sqlite3` to the CladeMers command.
This dependency is not needed when using `.dmp` taxonomy files.

## Unpack and compile

Transfer `CladeMers-0.13.0-source.zip` and its `.sha256` file to the HPC, then:

```bash
sha256sum -c CladeMers-0.13.0-source.zip.sha256
unzip CladeMers-0.13.0-source.zip
cd CladeMers-0.13.0
sha256sum -c SHA256SUMS

uname -m
g++ --version
blastdbcmd -version   # only needed for BLAST sequence input

make clean
make CXX=g++
./clademers --version
```

Clang can be used with `make CXX=clang++` instead. To build without Make:

```bash
c++ -O3 -std=c++17 -pthread clademers.cpp -lz -o clademers
```

`make clean` removes only `clademers`, `clademers.o`, and `clademers.dSYM`.
It succeeds when those artifacts are absent, and preserves source files,
databases, and analysis results. Run it before rebuilding for a different
architecture, compiler, or flags. Keep cleaning and building as separate commands;
do not run `make -j clean all`, whose goals can overlap.

Optional developer tests:

```bash
make test
# Include the genuine BLAST+ integration test; edit this directory:
make test-blast BLAST_BIN=/path/to/ncbi-blast/bin
# Also test against your existing BLAST taxonomy files (requires sqlite3):
make test-blast BLAST_BIN=/path/to/ncbi-blast/bin BLAST_TAXDB=/path/to/blast-taxonomy
```

`make test` skips the genuine BLAST+ integration test unless `BLAST_BIN` is set.
See `README.md` for a tiny synthetic database example and `VALIDATION.md` for
test coverage and platform limits.

## Extract Placozoa markers from NT

Run in a scheduler allocation with enough memory for all distinct target k-mers.
With 30 scanning workers, allocate 32 CPUs to leave room for the producer and
BLAST database reader. Memory and elapsed time for a full NT scan are not yet
measured. Replace the paths below with your actual taxonomy and full NT prefix:

If you have BLAST taxonomy files (`taxonomy4blast.sqlite3`, `taxdb.bti`, and
`taxdb.btd`), point to the directory containing all three:

```bash
./clademers \
  -blasttax /path/to/blast-taxonomy \
  -taxid 10226 -k 42 -t 30 \
  -blastdb /path/to/nt \
  -out /scratch/yourname/placozoa_k42
```

No separate `nodes.dmp`/`names.dmp` download is needed in this mode. BLAST's name
index sometimes omits higher-level taxa; CladeMers then displays the numeric
taxid, while still using the full SQLite tree to select the clade. Keep all
three taxonomy files unchanged during the run and use a standalone snapshot
without an active SQLite journal/WAL.

If you prefer NCBI taxonomy dump files, use the original interface instead:

```bash
./clademers \
  -taxnodes /path/to/taxonomy/nodes.dmp.gz \
  -taxnames /path/to/taxonomy/names.dmp.gz \
  -taxid 10226 \
  -k 42 \
  -t 30 \
  -blastdb /path/to/nt \
  -out /scratch/yourname/placozoa_k42
```

Plain `.dmp` files also work, and either dump can be compressed independently.
Compressed dumps require `gzip` on `PATH` or `-gzip /path/to/gzip`, including with
BLAST sequence input. Pass individual dump files, not the `taxdump.tar.gz`
archive itself. No decompressed copy is created.

Taxid **10226 selects Placozoa plus its descendants**. Use the full, real
NCBI taxonomy files matching your NT snapshot when possible; the bundled
`examples/nodes.dmp` and `examples/names.dmp` are tiny test fixtures.
Choose either `-blasttax` or the two `.dmp` options; do not combine them.
The `.dmp` files are distributed in
[`taxdump.tar.gz`](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz),
separately from BLAST's `taxdb.tar.gz`.
Use a complete, taxonomy-annotated version-5 nucleotide NT database prefix,
not a single numbered volume or a subset such as `core_nt`. Keep the database
unchanged during the run. The output directory must not already exist; its
parent directory must exist.

Every run now saves each marker's LCA taxid in its FASTA header and its taxid
and name in `markers.tsv`, plus `target_taxonomy.tsv`. Add **`-rd`** (equivalently
`--report-descendants`) to also write `descendant_summary.tsv` immediately.
You can instead create the report later without another NT scan:

```bash
./clademers -rd -results /scratch/yourname/placozoa_k42 \
  -out /scratch/yourname/placozoa_descendant_report
```

The later job needs only CladeMers and the completed run's `markers.fasta`,
`target_taxonomy.tsv`, `manifest.tsv`, and `SUCCESS`. No source reference or
external reader executable is needed; those four files can be moved off the
HPC. The summary includes every target-subtree taxon, including zero counts,
using the saved taxonomy snapshot. See README for the two count columns and
format validation. Report directories must be new. Older outputs cannot gain
missing all-occurrence LCA data from their representative source alone.

To also write one FASTA for each taxon with directly assigned markers, use
`-wsf` / `--write-separate-fastas`. It works independently of `-rd`:

```bash
./clademers -rd -wsf -results /scratch/yourname/placozoa_k42 \
  -out /scratch/yourname/placozoa_report_with_fastas
```

Files go under `taxon_fastas/taxid_<id>.fasta`; `taxon_fastas.tsv` maps taxids to
scientific names, direct counts and filenames. Each marker is written once,
under its exact LCA, without duplication into ancestors. Zero-count taxa have
no FASTA. Allow space for another copy of the combined marker FASTA and many
small files; the writer bounds open file handles. Saved 0.11.0 extractions can
be exported this way directly, without NT access. Add the same flag to an
initial extraction to write separate FASTAs as soon as it completes.

For compressed separate files, replace `-wsf` with `-wsf-gz` (long name:
`--write-separate-fastas-gz`):

```bash
./clademers -rd -wsf-gz -results /scratch/yourname/placozoa_k42 \
  -out /scratch/yourname/placozoa_report_with_gzip_fastas
```

This writes `taxon_fastas/taxid_<id>.fasta.gz` and records those paths in the
index. It implies `-wsf` and uses the already-linked zlib library at level 6;
no external gzip program is needed for output. The combined marker catalog
remains plain. Files may contain concatenated gzip members when the bounded
file cache reopens them; `gzip -dc` reads all members. Tiny files or frequent
reopening can reduce compression efficiency. The same flag works in an initial
extraction and with `-resume`.

Checkpoints save automatically every 300 seconds at safe record boundaries and
at phase transitions. Keep the result directory on shared persistent storage.
After a time limit or interruption, compile/install **0.13.0** for the compute
node architecture and run:

```bash
./clademers -resume /scratch/yourname/placozoa_k42 -t 30
```

Original input paths, taxonomy, taxid, k, and reader paths are restored. Do not
repeat those input arguments or `-out`; reader paths can be overridden when needed.
Use the same references and BLAST reader version. Both target loading and
background scanning resume directly by default. NT resumes at the saved OID,
revalidating only its saved defline rows; plain FASTA seeks to its byte offset;
gzip FASTA restores an indexed DEFLATE boundary and decompresses locally to the
saved record. No separate gzip indexing pass or decompressed copy is needed.
Taxonomy, maps, and saved candidate tables still need loading. Partial marker
output is rewritten from the completed scan checkpoint.

Keep reference files at their original paths, unchanged across submissions.
Direct resume compares paths, sizes, and modification/change timestamps for
FASTA or the BLAST alias/volume/filter files. Touching, replacing, or recopying
them can invalidate a checkpoint even if their content is identical. Use immutable
reference snapshots; these fast metadata checks are not full database checksums.
Gzip local replay depends on DEFLATE block size (normally short, but not a fixed
maximum). See README for details and the optional `-resume-mode replay` fallback.
The resume mode is chosen at initial run creation; LCA assignments are always
checkpointed. You may add `-rd`, `-wsf`, or `-wsf-gz` when resuming an incomplete
0.13.0 run. The compression choice is restored even if its flag is not repeated;
adding `-wsf-gz` enables compression for a run that began with plain output.
**Older checkpoints must be resumed with their original release.** Report-only
mode accepts completed 0.11.0 and newer extractions, rather than checkpoints.
The 0.13.0 checkpoint state format is 5 and saves the compression choice.

For jobs expected to span allocations, edit and submit
**`examples/nt_checkpoint.slurm`**. It uses a constant result path, automatically
chooses initial run or resume, and requests an advance signal:

```bash
sbatch examples/nt_checkpoint.slurm
# After it pauses or is killed, resubmit the same edited script:
sbatch examples/nt_checkpoint.slurm
```

The script requests `--signal=B:USR1@600`, forwards it from the batch shell to
CladeMers, and waits for saving to finish. The `B:` behavior is documented in
[Slurm's sbatch reference](https://slurm.schedmd.com/sbatch.html#OPT_signal).
Adjust the lead time to exceed your observed record/checkpoint-write latency;
large candidate tables can need more than ten minutes. A graceful pause exits
**75**, so Slurm may label it failed even though the checkpoint is valid. A hard
kill resumes from the last successful periodic save instead. The example does
not automatically resubmit jobs. Do not submit two jobs for the same result path
at once; CladeMers also locks the directory against concurrent writers.

`-checkpoint-seconds 600` changes the interval; `0` disables new saves. Saving
requires disk space: about 56 bytes per candidate plus source strings, with
room for two table copies during replacement/transition. After target loading,
the table is immutable and scan checkpoints contain compact flags/metadata.
Preserve both `checkpoint/state.bin` and `checkpoint/candidates.bin`, all inputs,
and the result directory. `make clean` preserves them. Older jobs started without
checkpoint support cannot be resumed retroactively. See README for restart and
filesystem requirements.

For Slurm, edit the executable, taxonomy, database, and result paths in
`examples/nt_blast_taxonomy.slurm` (BLAST taxonomy mode) or `examples/nt.slurm`
(`.dmp` mode), plus your site's account/partition and compiler/BLAST/SQLite
modules as needed, then submit the selected script:

```bash
sbatch examples/nt_blast_taxonomy.slurm
```

## Use FASTA instead of a BLAST database

FASTA input does not need BLAST+. Use a plain or gzip-compressed regular FASTA file with
both target and background references. If headers contain explicit `taxid=`,
`taxids=`, or Kraken taxid annotations:

```bash
./clademers -fasta /path/to/references.fasta.gz \
  -blasttax /path/to/blast-taxonomy \
  -taxid 10226 -k 42 -t 30 -out /scratch/yourname/placozoa_fasta_k42
```

If taxids are in a separate two-column mapping file, add
`-taxidmap /path/to/sequence_taxids.tsv` or
`-taxidmap /path/to/sequence_taxids.tsv.gz`. The key must match the exact first
whitespace-delimited FASTA ID. See `README.md` for supported header/map formats,
handling of multiple/missing taxids, and conflicts. NCBI's four-column
accession-to-taxid tables must be converted to this two-column format first.

The taxonomy options are independent of sequence input. Replace `-blasttax`
with `-taxnodes /path/to/nodes.dmp.gz -taxnames /path/to/names.dmp.gz` to use dumps.
With all text inputs plain, FASTA plus dumps needs no external runtime executable.
Plain or gzip FASTA with plain maps plus BLAST taxonomy needs only the `sqlite3`
CLI as an external program. Compressed maps/dumps and replay-mode gzip FASTA
also need `gzip` on `PATH`, or `-gzip /path/to/gzip`. Keep `-blastdbcmd` off FASTA commands.

The same command accepts `.fasta`, `.fa`, `.fasta.gz`, and `.fa.gz`; gzip is
detected from file contents. No decompressed copy is written to disk. Indexed
gzip FASTA uses the zlib library linked during compilation. FASTA,
mapping files, nodes.dmp, and names.dmp can all be compressed independently;
the decompressed map is still
two-column text. Decompression errors
invalidate the run rather than accepting a partial reference scan. For example,
with every text input compressed:

```bash
./clademers -fasta /path/to/references.fasta.gz \
  -taxidmap /path/to/sequence_taxids.tsv.gz \
  -taxnodes /path/to/nodes.dmp.gz -taxnames /path/to/names.dmp.gz \
  -taxid 10226 -k 42 -t 30 -out /scratch/yourname/placozoa_k42
```

FASTA is read twice, decompressing on each pass when needed; one complete
sequence and the optional mapping file are buffered in RAM in addition to
target k-mers. The map is streamed again for fingerprint verification, including
decompression when needed, without allocating a second map. Taxonomy dumps are
also streamed again for verification without rebuilding the tree. Gzip compression
saves disk space; loaded map/tree RAM requirements remain the same.
Decompressed content fingerprints and all reader exit statuses are checked
before success. Keep FASTA, map, and taxonomy files unchanged during the run.
Both `-blastdb` and `-fasta` retain the same exact k-mer filtering rules.

The Slurm examples' 64 GB memory and two-day time requests are placeholders, not
measured requirements. Lower `-t` if your allocation has fewer CPUs.

Progress percentages, elapsed time, and approximate ETA go to stderr and appear
automatically in `clademers-<jobid>.log` with these Slurm examples. The default
update interval is 30 seconds; add `-progress-seconds 10` to change it. Watch the
file with `tail -f clademers-<jobid>.log` (replace `<jobid>` with the actual ID).
Updates occur at processing boundaries, so long records/read stalls can delay them.

Percentages and ETAs apply to the named operation, not the entire job. The NT
background scan uses BLAST's reported base count; the second FASTA pass also
has a known content size when compressed. First gzip passes and BLAST target
extraction report `ETA unavailable` until their totals are known at completion.
Snapshot checks and marker output follow the sequence scan. See README for the
complete progress accounting and unknown-total behavior.

Accept a result only when `SUCCESS` exists and `manifest.tsv` has
`status` equal to `complete`. The full marker catalog is in `markers.fasta`
and `markers.tsv`. Specificity means present in target references and absent
from screened background references in this database snapshot; it does not
establish universal biological specificity or organism presence in a sample.
