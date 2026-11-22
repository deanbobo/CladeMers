# Publishing CladeMers

Repository: <https://github.com/deanbobo/CladeMers>

The initial repository includes the C++ source, Makefile, documentation, MIT
license, citation metadata, synthetic examples, tests, and GitHub Actions
configuration. Local research drafts, session notes, reference data, compiled
binaries, and historical ZIPs are excluded by `.gitignore`.

## First push

If starting from an unpacked source ZIP in a new directory, initialize Git and
stage the source files first (skip this block for the already-prepared checkout):

```bash
git init -b main
git remote add origin https://github.com/deanbobo/CladeMers.git
git add .gitattributes .gitignore .github CITATION.cff LICENSE \
  CODE_WALKTHROUGH.md CONTRIBUTING.md GITHUB_SETUP.md HPC.md \
  Makefile README.md VALIDATION.md clademers.cpp examples tests
```

From the prepared project directory, inspect the staged files and confirm your
Git author identity:

```bash
git status
git diff --cached --stat
git config user.name
git config user.email
```

If you need to set an identity for this repository, use `git config user.name`
and `git config user.email` with your chosen values. A GitHub-provided noreply
email can be used for commit privacy.

The local repository uses `main` and the HTTPS origin above. Create the initial
commit and push it when ready:

```bash
git commit -m "Prepare CladeMers 0.13.0 for public release"
git push -u origin main
```

Authenticate using your configured Git credential manager or GitHub's supported
Git authentication method. Do not put access tokens in repository files or
remote URLs. See GitHub's
[instructions for importing local code](https://docs.github.com/en/migrations/importing-source-code/using-the-command-line-to-import-source-code/adding-locally-hosted-code-to-github).

## After the first push

1. Open the **Actions** tab and check both Linux AMD64 and ARM64 test jobs.
   The workflow has been prepared locally; hosted success requires an actual run.
2. Set the repository description to: “Exact lineage-specific k-mer discovery
   from BLAST databases or FASTA, with taxonomy-aware descendant reports.”
   Suggested topics: `bioinformatics`, `kmers`, `taxonomy`, `edna`, `fasta`, `blast`.
3. Check the README, MIT license, and **Cite this repository** entry on GitHub.
4. When the checks pass, create a GitHub Release tagged `v0.13.0`. Describe the
   supported inputs, saved LCA assignments, gzip taxon exports, and checkpoint
   compatibility. State that full NT/HPC performance validation is still pending.

GitHub supplies source archives for a tag. Optional source ZIPs and their
checksums can be attached as release assets; keep them out of Git history.
An archived release DOI can be added to `CITATION.cff` later after one is issued.

The tests create small local databases; they do not download NT. Optional
taxonomy-snapshot, prior-release, and cross-architecture checkpoint fixtures
are described in [CONTRIBUTING.md](CONTRIBUTING.md).
