# Descendant report preview

This is a synthetic example, not a survey of Placozoa or an NT marker catalog.
The `Example ...` taxa and their IDs are fictional. All sequences are synthetic.
The six distinct canonical 42-mers demonstrate shared and private assignments.

From the source-package root, run an extraction **without** `-rd`:

```bash
./clademers -fasta examples/descendant_demo/references.fasta \
  -taxnodes examples/descendant_demo/nodes.dmp \
  -taxnames examples/descendant_demo/names.dmp \
  -taxid 10226 -k 42 -out descendant_demo_results
```

Then create the report from those saved results:

```bash
./clademers -rd -wsf -results descendant_demo_results -out descendant_demo_report
```

The bundled `markers.fasta`, `markers.tsv`, and `descendant_summary.tsv` were
generated with these commands and are expected outputs. The accompanying
`taxon_fastas.tsv` and `taxon_fastas/` show the separate-FASTA export. The unit tests verify
them. Only the generated result directory includes the manifest and saved
taxonomy required by `-results`; the loose preview files are examples.

| Clade | `assigned_kmers` | `assigned_kmers_incl_descendants` |
| --- | ---: | ---: |
| Placozoa | 1 | 6 |
| Example lineage A | 1 | 4 |
| Example species A1 | 2 | 2 |
| Example species A2 | 1 | 1 |
| Example lineage B | 1 | 1 |
| Example unsampled lineage | 0 | 0 |

`assigned_kmers` counts assignments directly to the node.
`assigned_kmers_incl_descendants` counts assignments to the node plus all its
descendants.

One marker shared by A1 and A2 is assigned to lineage A. Another shared by A1
and lineage B is assigned to Placozoa. The remaining four markers are assigned
to A1 (two), A2 (one), and B (one). Markers appear once in the catalog; subtree
counts include markers assigned farther down the tree and therefore overlap.


The separate FASTAs are:

```text
taxon_fastas/taxid_10226.fasta    # 1 marker assigned directly to Placozoa
taxon_fastas/taxid_9000001.fasta  # 1 marker shared by A1 and A2
taxon_fastas/taxid_9000002.fasta  # 2 markers assigned directly to A1
taxon_fastas/taxid_9000003.fasta  # 1 marker assigned directly to A2
taxon_fastas/taxid_9000004.fasta  # 1 marker assigned directly to B
```

The zero-count lineage has no file. Child markers are not repeated in ancestor
files. The index records names and counts; FASTA headers retain the original
marker IDs and `taxid=<LCA>` descriptions. `-wsf` can also be used without `-rd`.

To write the same per-taxon content as gzip-compressed FASTAs:

```bash
./clademers -rd -wsf-gz -results descendant_demo_results -out descendant_demo_gzip
gzip -dc descendant_demo_gzip/taxon_fastas/taxid_9000002.fasta.gz
```

`-wsf-gz` (or `--write-separate-fastas-gz`) implies `-wsf`. The index uses
`.fasta.gz` paths, and decompressing each file reproduces the corresponding
bundled plain FASTA. The writer uses zlib directly; the second command merely
illustrates how to view one output. No NT scan or external gzip process is
needed to generate the compressed export. These tiny demo files may be larger
after compression.
