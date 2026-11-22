# Contributing to CladeMers

Bug reports, reproducible examples, documentation improvements, and code changes
are welcome. Please discuss substantial changes in a GitHub issue first.

## Report a problem

Include `./clademers --version`, your operating system and CPU architecture,
compiler version, exact command, and the relevant error output. State whether
the run uses FASTA or BLAST and `.dmp` or BLAST taxonomy. For BLAST problems,
include the BLAST+ version. For restart problems, include the original command,
resume command, and checkpoint phase.

A small synthetic input with expected and actual output is especially useful.
Please remove private paths and unpublished sequences before posting logs or
examples. There is no need to attach NT or an entire reference database.

## Build and test

Use a C++17 compiler and zlib development headers/library:

```bash
make clean
make
make test
```

Python 3's standard library runs the tests; it is not required to compile or run
CladeMers. `make test` uses temporary synthetic datasets and does not download
NT. Install the `gzip` and `sqlite3` command-line tools to enable those test
cases. To include integration with BLAST+ 2.15 or newer:

```bash
make test-blast BLAST_BIN=/path/to/ncbi-blast/bin
```

Optional fixtures enable additional cases:

- `BLAST_TAXDB`: an extracted NCBI BLAST taxonomy snapshot.
- `CLADEMERS_OLD_BINARY`: a CladeMers 0.11.0 executable, for portable-result compatibility.
- `CLADEMERS_RESUME_BINARY`: a runnable second-architecture build of the current
  release, for checkpoint exchange. It must be executable on the test host.
- `CLADEMERS_BINARY`: an alternative build to test, such as a sanitizer build.

Missing optional tools or fixtures produce explicit test skips. See
[VALIDATION.md](VALIDATION.md) for recorded validation and its limits.

## Submit a change

Keep changes focused and explain the user-visible behavior. Add or update tests
when changing marker assignment, input handling, saved formats, compression, or
checkpoint recovery. Prefer small inputs with independently calculated answers.
Update the help text and README for changes to command-line options or outputs.

Preserve canonical reverse-complement handling, conservative treatment of unknown
or mixed taxonomy, and the distinction between a direct LCA assignment and counts
including descendants. Changes to serialized state must explicitly address
checkpoint versions and portable-result compatibility.

Open a pull request with the problem, the resulting behavior, and the checks you
ran. Include any relevant test skips. Avoid committing executables, reference
databases, result directories, or release ZIPs.

## Automated checks

The GitHub Actions workflow builds with GCC and runs the test suite on Linux
AMD64 and ARM64. It installs BLAST+ 2.17.0 for genuine database integration;
the large NCBI taxonomy snapshot, prior-release executable, and cross-architecture
checkpoint fixture remain optional local checks. The workflow does not download
NT. The first hosted run is needed to establish CI results for a new repository.
Runner labels follow the [GitHub-hosted runner reference](https://docs.github.com/en/actions/reference/runners/github-hosted-runners).

Maintainers can follow [GITHUB_SETUP.md](GITHUB_SETUP.md) for the first push and
release steps.
