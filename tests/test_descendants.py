#!/usr/bin/env python3
"""Independent string/set oracle for most-specific descendant assignments."""
import csv
import gzip
import os
from pathlib import Path
import random
import resource
import shutil
import signal
import subprocess
import unittest
import zlib

import test_clademers as core
import test_checkpoint as restarts


class DescendantTests(unittest.TestCase):
    tearDown = restarts.CheckpointTests.tearDown
    configure = restarts.CheckpointTests.configure
    run_ok = restarts.CheckpointTests.run_ok
    run_bad = restarts.CheckpointTests.run_bad
    launch = restarts.CheckpointTests.launch
    wait_log = restarts.CheckpointTests.wait_log
    interrupt = restarts.CheckpointTests.interrupt
    compare = restarts.CheckpointTests.compare
    manifest = restarts.CheckpointTests.manifest

    def setUp(self):
        restarts.CheckpointTests.setUp(self)
        self.parents = {1: 1, 33208: 1, 10226: 33208, 10227: 10226,
                        10228: 10226, 9000001: 10227, 9000002: 10227,
                        9000003: 10228, 9606: 33208, 2: 1}
        (self.root/'nodes.dmp').write_text(''.join(
            f'{node}\t|\t{parent}\t|\tno rank\t|\n' for node, parent in self.parents.items()))
        (self.root/'names.dmp').write_text(''.join(
            f'{node}\t|\t{name}\t|\t\t|\tscientific name\t|\n'
            for node, name in ((n, 'taxon '+str(n)) for n in self.parents)))
        rng = random.Random(173)
        a, b, c, genus, root, broad, outside, unknown, multi, mixed = [
            ''.join(rng.choices('ACGT', k=90)) for _ in range(10)]
        self.entries = [
            ('a', 'N'.join((a, genus, root, broad, outside, unknown)), '9000001', 'N/A'),
            ('b', 'N'.join((b, core.reverse_complement(genus), root)), '9000002', 'N/A'),
            ('c', c+'N'+root, '10228', 'N/A'),
            ('broad', broad, '10227', 'N/A'),
            ('reverse', core.reverse_complement(a).lower(), '9000001', 'N/A'),
            ('multi', multi, '9000001;9000002', 'N/A'),
            ('mixed', mixed, '9000001;9606', 'N/A'),
            ('unknown_match', unknown, 'N/A', 'N/A'),
            ('outside_match', outside, '9606', 'N/A'),
            ('zero', 'A'*90, '0', 'N/A')]

    def lineage(self, taxid):
        result = []
        while taxid in self.parents:
            result.append(taxid)
            if self.parents[taxid] == taxid:
                break
            taxid = self.parents[taxid]
        return result

    def annotations(self, entry):
        fields = [entry[2]]
        if entry[3] not in ('N/A', '0', ''):
            fields.append(entry[3])
        return {int(v) if v not in ('N/A', '') else 0
                for field in fields for v in field.replace(',', ';').split(';')}

    def oracle(self, entries, k, target=10226):
        observations = {}
        for entry in entries:
            for key in core.kmers(entry[1], k):
                observations.setdefault(key, set()).update(self.annotations(entry))
        result = {}
        for key, ids in observations.items():
            if all(target in self.lineage(t) for t in ids):
                paths = [self.lineage(t) for t in ids]
                common = set(paths[0]).intersection(*map(set, paths[1:]))
                result[key] = next(t for t in paths[0] if t in common)
        return result

    def common(self, output, target=10226, k=42, threads=1):
        return ['-taxnodes', str(self.root/'nodes.dmp'), '-taxnames', str(self.root/'names.dmp'),
                '-taxid', str(target), '-k', str(k), '-t', str(threads), '-out', str(output),
                '-checkpoint-seconds', '1']

    def blast_args(self, output, entries=None, target=10226, k=42, flag='-rd'):
        entries = self.entries if entries is None else entries
        inside = [e for e in entries if any(target in self.lineage(t) for t in self.annotations(e))]
        (self.db/'target.tsv').write_text(''.join(core.row(*e) for e in inside))
        (self.db/'all.tsv').write_text(''.join(core.row(*e) for e in entries))
        self.ids = [str(t) for t in self.parents if target in self.lineage(t)]
        self.configure(expected_taxids=self.ids)
        return self.common(output, target, k)+['-blastdb', str(self.db), '-blastdbcmd', str(self.fake),
                '-resume-mode', 'replay']+([flag] if flag else [])

    def table(self, output, name='markers.tsv'):
        with (output/name).open() as handle:
            return list(csv.DictReader(handle, delimiter='\t'))

    def check_report(self, output, entries=None, k=42, target=10226):
        expected = self.oracle(self.entries if entries is None else entries, k, target)
        markers = self.table(output)
        self.assertEqual({r['sequence']: int(r['specificity_taxid']) for r in markers}, expected)
        self.assertEqual(len(markers), len(expected))
        summary = self.table(output, 'descendant_summary.tsv')
        self.assertEqual([int(r['taxid']) for r in summary],
                         sorted(t for t in self.parents if target in self.lineage(t)))
        for row in summary:
            taxid = int(row['taxid'])
            self.assertEqual(int(row['parent_taxid']), self.parents[taxid])
            self.assertEqual(int(row['assigned_kmers']), sum(t == taxid for t in expected.values()))
            self.assertEqual(int(row['assigned_kmers_incl_descendants']),
                             sum(taxid in self.lineage(t) for t in expected.values()))
        self.assertEqual(self.manifest(output)['report_descendants'], 'yes')
        self.assertTrue((output/'SUCCESS').exists())
        return expected

    def test_lca_always_saved_aliases_and_rejected_overrides(self):
        baseline = self.root/'baseline'
        self.run_ok(self.blast_args(baseline, flag=None))
        self.assertEqual({r['sequence']: int(r['specificity_taxid']) for r in self.table(baseline)}, self.oracle(self.entries, 42))
        headers = (baseline/'markers.fasta').read_text().splitlines()[::2]
        self.assertEqual(headers, ['>'+r['marker_id']+' taxid='+r['specificity_taxid'] for r in self.table(baseline)])
        self.assertTrue((baseline/'target_taxonomy.tsv').exists())
        self.assertFalse((baseline/'descendant_summary.tsv').exists())
        for i, flag in enumerate(('-rd', '--report-descendants')):
            output = self.root/('alias'+str(i))
            self.run_ok(self.blast_args(output, flag=flag))
            self.check_report(output)
            self.assertEqual((output/'markers.fasta').read_bytes(), (baseline/'markers.fasta').read_bytes())
            self.assertEqual((output/'markers.tsv').read_bytes(), (baseline/'markers.tsv').read_bytes())
        args = self.blast_args(self.root/'invalid')
        self.run_bad(args+['--report-descendants'], 'Repeated option')
        self.run_bad(['-resume', str(baseline), '-rd'], 'already complete')

    def test_oracle_bit_boundaries_threads_and_blast_leaf_annotations(self):
        entries = self.entries + [
            ('leaves', self.entries[0][1], '9000001', '9000002'),
            ('broad_leaf', self.entries[2][1], '10226', '10228'),
            ('mixed_leaf', self.entries[5][1], '9000001', '9606'),
            ('unknown_leaf', self.entries[3][1], '9000001', '9999999')]
        for k in (1, 7, 32, 33, 42, 64):
            previous = None
            for threads in (1, 4):
                with self.subTest(k=k, threads=threads):
                    output = self.root/f'k{k}-t{threads}'
                    args = self.blast_args(output, entries, k=k)
                    args[args.index('-t')+1] = str(threads)
                    self.run_ok(args)
                    self.check_report(output, entries, k)
                    if previous:
                        self.assertEqual((output/'markers.tsv').read_bytes(), (previous/'markers.tsv').read_bytes())
                    previous = output

    def test_assignments_equal_separate_descendant_runs(self):
        output = self.root/'root'
        self.run_ok(self.blast_args(output))
        assignments = self.check_report(output)
        self.assertEqual(set(assignments.values()), {10226, 10227, 10228, 9000001, 9000002})
        for target in (10227, 10228, 9000001, 9000002):
            child = self.root/str(target)
            self.run_ok(self.blast_args(child, target=target, flag=None))
            self.assertEqual({r['sequence'] for r in self.table(child)},
                {key for key, taxid in assignments.items() if target in self.lineage(taxid)})

    @unittest.skipUnless(core.GZIP_BINARY and core.SQLITE3_BINARY, 'gzip and sqlite3 required')
    def test_fasta_compressed_text_inputs_both_taxonomies_and_missing_names(self):
        taxonomy = core.write_blast_taxonomy(self.root/'blast-taxonomy', self.root)
        for key in ('nodes', 'names'):
            (self.root/(key+'.dmp.gz')).write_bytes(gzip.compress((self.root/(key+'.dmp')).read_bytes()))
        for mapped, compressed, blasttax in ((False, False, False), (True, True, False), (False, True, True)):
            output = self.root/f'fasta-{mapped}-{compressed}-{blasttax}'
            content = ''.join('>'+e[0]+('' if mapped else ' taxids='+e[2])+'\n'+e[1]+'\n' for e in self.entries)
            fasta = self.root/(output.name+'.fa'+('.gz' if compressed else ''))
            fasta.write_bytes(gzip.compress(content.encode()) if compressed else content.encode())
            args = self.common(output)+['-fasta', str(fasta), '-rd']
            if mapped:
                mapping = self.root/'mapping.tsv.gz'
                mapping.write_bytes(gzip.compress(''.join(e[0]+'\t'+e[2]+'\n' for e in self.entries).encode()))
                args += ['-taxidmap', str(mapping)]
            if blasttax:
                args = args[4:]+['-blasttax', str(taxonomy), '-sqlite3', core.SQLITE3_BINARY]
            elif compressed:
                args[1] += '.gz'; args[3] += '.gz'
            args += ['-gzip', core.GZIP_BINARY]
            self.run_ok(args)
            self.check_report(output)
            self.assertTrue(all(r['specificity_name'] == 'taxon '+r['specificity_taxid'] for r in self.table(output)))
        # BLAST's name database can omit an internal node without losing membership.
        (self.root/'names.dmp').write_text(''.join(line for line in (self.root/'names.dmp').read_text().splitlines(True)
                                                  if not line.startswith('10227\t')))
        missing = core.write_blast_taxonomy(self.root/'missing-names', self.root)
        output = self.root/'fallback'
        args = self.blast_args(output)[4:]+['-blasttax', str(missing), '-sqlite3', core.SQLITE3_BINARY]
        self.run_ok(args)
        self.check_report(output)
        self.assertTrue(all(r['specificity_name'] == 'taxid 10227' for r in self.table(output)
                            if r['specificity_taxid'] == '10227'))

    def test_empty_summary_leaf_target_and_root_target(self):
        entries = [('a', 'ACGT'*30, '9000001', 'N/A'), ('b', 'ACGT'*30, '9606', 'N/A')]
        output = self.root/'empty'
        self.run_ok(self.blast_args(output, entries))
        self.assertEqual(self.check_report(output, entries), {})
        for target in (1, 9000001):
            output = self.root/('target'+str(target))
            self.run_ok(self.blast_args(output, self.entries, target=target))
            self.check_report(output, target=target)

    def test_replay_checkpoints_target_background_and_output(self):
        baseline = self.root/'baseline'
        self.run_ok(self.blast_args(baseline))
        for phase in ('target loading', 'background scan', 'output'):
            output = self.root/phase.replace(' ', '-')
            if phase == 'output':
                shutil.copytree(baseline, output)
                (output/'SUCCESS').unlink()
                (output/'descendant_summary.tsv').write_text('partial garbage')
            else:
                args = self.blast_args(output, flag=None)
                self.configure(expected_taxids=self.ids, **{
                    'target_row_delay' if phase == 'target loading' else 'row_delay': 0.3})
                process, log = self.launch(args)
                self.wait_log(process, log, 'Saved '+phase+r' after [1-9][0-9]* records', regex=True)
                self.interrupt(process, log, signal.SIGKILL)
                self.configure(expected_taxids=self.ids)
            self.run_ok(['-resume', str(output), '-t', '4', '-rd'])
            self.compare(baseline, output)
            self.check_report(output)
            self.assertEqual((baseline/'descendant_summary.tsv').read_bytes(), (output/'descendant_summary.tsv').read_bytes())

    def test_direct_plain_and_indexed_gzip_target_resume(self):
        rng = random.Random(324)
        motif = ''.join(rng.choices('ACGT', k=256))
        oracle_entries = [('a', motif*2, '9000001', 'N/A'), ('b', motif*2, '9000002', 'N/A'),
                          ('background', 'C'*42, '9606', 'N/A')]
        other = Path(os.environ.get('CLADEMERS_RESUME_BINARY', str(core.BINARY)))
        for compressed in (False, True):
            path = self.root/('large.fa.gz' if compressed else 'large.fa')
            with (gzip.open(path, 'wt', compresslevel=1) if compressed else path.open('w')) as handle:
                for i in range(128):
                    handle.write(f'>s{i} taxid={9000002 if i == 127 else 9000001}\n'+motif*4096+'\n')
                handle.write('>background taxid=9606\n'+'C'*42+'\n')
            baseline = self.root/('direct-baseline'+str(compressed))
            self.run_ok(self.common(baseline)+['-fasta', str(path), '-rd'])
            output = self.root/('direct-paused'+str(compressed))
            process, log = self.launch(self.common(output)+['-fasta', str(path)])
            self.wait_log(process, log, r'Saved target loading after [1-9][0-9]* records', regex=True)
            self.interrupt(process, log, signal.SIGKILL)
            result = self.run_ok(['-resume', str(output), '-t', '4', '-rd'], binary=other)
            self.assertIn('Indexed gzip seek' if compressed else 'Plain FASTA seek', result.stderr)
            self.compare(baseline, output)
            self.check_report(output, oracle_entries)
            self.assertEqual((baseline/'descendant_summary.tsv').read_bytes(), (output/'descendant_summary.tsv').read_bytes())

    @unittest.skipUnless(os.environ.get('BLAST_BIN'), 'real BLAST+ required')
    def test_real_blast_matches_fasta_all_occurrences(self):
        entries = [e for e in self.entries if len(self.annotations(e)) == 1 and 0 not in self.annotations(e)]
        # Equal sequences assigned to different taxa must contribute both labels.
        entries += [('same_sequence_other_species', entries[0][1], '9000002', 'N/A')]
        fasta = self.root/'real.fa'; mapping = self.root/'real.tsv'
        fasta.write_text(''.join('>'+e[0]+'\n'+e[1]+'\n' for e in entries))
        mapping.write_text(''.join(e[0]+'\t'+e[2]+'\n' for e in entries))
        database = self.root/'real-db'; bindir = Path(os.environ['BLAST_BIN'])
        made = subprocess.run([str(bindir/'makeblastdb'), '-in', str(fasta), '-dbtype', 'nucl',
            '-parse_seqids', '-taxid_map', str(mapping), '-out', str(database)], capture_output=True, text=True)
        self.assertEqual(made.returncode, 0, made.stderr)
        for mode in ('blast', 'fasta'):
            output = self.root/('real-'+mode)
            args = (['-blastdb', str(database), '-blastdbcmd', str(bindir/'blastdbcmd')] if mode == 'blast'
                    else ['-fasta', str(fasta), '-taxidmap', str(mapping)])
            self.run_ok(self.common(output)+args+['-rd'])
            self.check_report(output, entries)

    def saved_report(self, source, output=None, success=True, extra=None):
        existed = output is not None and output.exists()
        args = [str(core.BINARY), '-rd', '-results', str(source)]
        if output is not None:
            args += ['-out', str(output)]
        result = subprocess.run(args+(extra or []), capture_output=True, text=True, cwd=self.root,
                                env=dict(os.environ, PATH=str(self.root/'no-programs')), timeout=30)
        self.assertEqual(result.returncode == 0, success, result.stderr)
        if success:
            self.assertIn('no reference scan', result.stderr)
            report = output or source/'descendant_report'
            self.assertTrue((report/'SUCCESS').exists())
            self.assertEqual(self.manifest(report)['mode'], 'descendant_report')
            self.assertFalse((report/'commands.log').exists())
        elif output is not None and not existed:
            self.assertFalse((output/'SUCCESS').exists())
        return result

    def test_later_report_without_original_inputs_tsv_or_checkpoints(self):
        immediate = self.root/'immediate'
        self.run_ok(self.blast_args(immediate))
        baseline = self.root/'baseline'
        self.run_ok(self.blast_args(baseline, flag=None))
        expected = (immediate/'descendant_summary.tsv').read_bytes()
        moved = self.root/'moved portable results'; moved.mkdir()
        for name in ('markers.fasta', 'target_taxonomy.tsv', 'manifest.tsv', 'SUCCESS'):
            shutil.copyfile(baseline/name, moved/name)
        before = {p.name: p.read_bytes() for p in moved.iterdir()}
        shutil.rmtree(self.db)
        self.fake.unlink()
        (self.root/'nodes.dmp').unlink(); (self.root/'names.dmp').unlink()
        self.saved_report(moved)  # Default output subdirectory; no reader executables in PATH.
        self.assertEqual((moved/'descendant_report/descendant_summary.tsv').read_bytes(), expected)
        for name, data in before.items():
            self.assertEqual((moved/name).read_bytes(), data)
        output = self.root/'separate-report'
        self.saved_report(moved, output)
        self.assertEqual((output/'descendant_summary.tsv').read_bytes(), expected)
        self.saved_report(moved, output, success=False)  # Existing outputs cannot be overwritten.
        self.assertEqual((output/'descendant_summary.tsv').read_bytes(), expected)

    @unittest.skipUnless(core.GZIP_BINARY and core.SQLITE3_BINARY, 'gzip and sqlite3 required')
    def test_later_report_from_compressed_fasta_and_blast_taxonomy(self):
        taxonomy = core.write_blast_taxonomy(self.root/'blast-taxonomy', self.root)
        path = self.root/'references.fa.gz'
        path.write_bytes(gzip.compress(''.join('>'+e[0]+' taxids='+e[2]+'\n'+e[1]+'\n' for e in self.entries).encode()))
        baseline = self.root/'baseline'
        args = self.common(baseline)[4:]+['-blasttax', str(taxonomy), '-sqlite3', core.SQLITE3_BINARY, '-fasta', str(path)]
        self.run_ok(args)
        expected = self.oracle(self.entries, 42)
        self.assertEqual({r['sequence']: int(r['specificity_taxid']) for r in self.table(baseline)}, expected)
        path.unlink(); shutil.rmtree(taxonomy)
        self.saved_report(baseline)
        summary = self.table(baseline/'descendant_report', 'descendant_summary.tsv')
        self.assertEqual(sum(int(r['assigned_kmers']) for r in summary), len(expected))
        self.assertEqual({int(r['taxid']): int(r['assigned_kmers']) for r in summary},
                         {t: sum(v == t for v in expected.values()) for t in self.parents if 10226 in self.lineage(t)})

    def test_later_report_empty_catalog_and_root_or_leaf(self):
        for target in (1, 10226, 9000001):
            entries = [('a', 'ACGT'*30, '9000001', 'N/A'), ('b', 'ACGT'*30, '9606', 'N/A'),
                       ('unknown', 'ACGT'*30, 'N/A', 'N/A')]
            output = self.root/str(target)
            self.run_ok(self.blast_args(output, entries, target=target))
            self.assertEqual((output/'markers.fasta').read_bytes(), b'')
            self.saved_report(output)
            self.assertEqual((output/'descendant_summary.tsv').read_bytes(),
                             (output/'descendant_report/descendant_summary.tsv').read_bytes())

    def test_later_report_rejects_incomplete_legacy_and_damaged_results(self):
        baseline = self.root/'baseline'
        self.run_ok(self.blast_args(baseline, flag=None))
        cases = {
            'missing-success': lambda p: (p/'SUCCESS').unlink(),
            'legacy': lambda p: self.edit_manifest(p, 'results_format', None),
            'incomplete': lambda p: self.edit_manifest(p, 'status', 'running'),
            'count': lambda p: self.edit_manifest(p, 'retained_kmers', '0'),
            'missing-fasta': lambda p: (p/'markers.fasta').unlink(),
            'missing-tree': lambda p: (p/'target_taxonomy.tsv').unlink(),
            'changed-name': lambda p: (p/'target_taxonomy.tsv').write_bytes((p/'target_taxonomy.tsv').read_bytes().replace(b'taxon', b'alter', 1)),
            'truncated': lambda p: (p/'markers.fasta').write_bytes((p/'markers.fasta').read_bytes()[:-1]),
            'no-lca': lambda p: (p/'markers.fasta').write_bytes((p/'markers.fasta').read_bytes().replace(b' taxid=', b' label=', 1)),
            'outside-lca': lambda p: (p/'markers.fasta').write_bytes((p/'markers.fasta').read_bytes().replace(b'taxid=10227', b'taxid=9606', 1)),
            'changed-lca': lambda p: (p/'markers.fasta').write_bytes((p/'markers.fasta').read_bytes().replace(b'taxid=10227', b'taxid=10226', 1)),
        }
        for name, mutate in cases.items():
            with self.subTest(name=name):
                copy = self.root/name; shutil.copytree(baseline, copy)
                mutate(copy)
                self.saved_report(copy, self.root/('report-'+name), success=False)

    def edit_manifest(self, directory, key, value):
        fields = self.manifest(directory)
        if value is None:
            fields.pop(key, None)
        else:
            fields[key] = str(value)
        (directory/'manifest.tsv').write_text('key\tvalue\n'+''.join(k+'\t'+v+'\n' for k, v in fields.items()))

    def test_later_report_validates_saved_tree_and_cli(self):
        baseline = self.root/'baseline'
        self.run_ok(self.blast_args(baseline, flag=None))
        for name, old, new in (
                ('cycle', b'10227\t10226', b'10227\t9000001'),
                ('missing-parent', b'9000001\t10227', b'9000001\t9999999'),
                ('bad-root', b'10226\t33208', b'10226\t9000001'),
                ('duplicate', b'9000002\t10227', b'9000001\t10227')):
            with self.subTest(name=name):
                copy = self.root/name; shutil.copytree(baseline, copy)
                data = (copy/'target_taxonomy.tsv').read_bytes().replace(old, new)
                (copy/'target_taxonomy.tsv').write_bytes(data)
                self.edit_manifest(copy, 'target_taxonomy_fnv1a64', core.fnv1a64(data))
                self.edit_manifest(copy, 'target_taxonomy_bytes', len(data))
                self.saved_report(copy, self.root/('report-'+name), success=False)
        for extra in (['-fasta', 'anything'], ['-blastdb', 'nt'], ['-taxid', '10226'],
                      ['-taxnodes', 'nodes.dmp'], ['-resume', str(baseline)], ['-rd']):
            self.saved_report(baseline, self.root/'invalid-cli', success=False, extra=extra)
        self.run_bad(['-results', str(baseline)], 'requires -rd')

    def test_bundled_preview_outputs(self):
        demo = core.ROOT/'examples/descendant_demo'
        output = self.root/'demo'
        self.run_ok(['-fasta', str(demo/'references.fasta'), '-taxnodes', str(demo/'nodes.dmp'),
                     '-taxnames', str(demo/'names.dmp'), '-taxid', '10226', '-out', str(output)])
        self.saved_report(output, extra=['-wsf'])
        for name in ('markers.fasta', 'markers.tsv'):
            self.assertEqual((demo/name).read_bytes(), (output/name).read_bytes())
        self.assertEqual((demo/'descendant_summary.tsv').read_bytes(),
                         (output/'descendant_report/descendant_summary.tsv').read_bytes())
        self.assertEqual((demo/'taxon_fastas.tsv').read_bytes(),
                         (output/'descendant_report/taxon_fastas.tsv').read_bytes())
        for path in (demo/'taxon_fastas').glob('*.fasta'):
            self.assertEqual(path.read_bytes(), (output/'descendant_report/taxon_fastas'/path.name).read_bytes())

    def check_separate(self, output, catalog, compressed=False):
        lines = catalog.read_text().splitlines(keepends=True)
        groups = {}
        for header, sequence in zip(lines[::2], lines[1::2]):
            taxid = int(header.split('taxid=')[1])
            groups.setdefault(taxid, []).append(header+sequence)
        index = self.table(output, 'taxon_fastas.tsv')
        suffix = '.fasta.gz' if compressed else '.fasta'
        self.assertEqual([int(r['taxid']) for r in index], sorted(groups))
        names = {int(r['taxid']): r['scientific_name'] for r in self.table(catalog.parent, 'target_taxonomy.tsv')}
        for row in index:
            taxid = int(row['taxid'])
            self.assertEqual(row['scientific_name'], names[taxid])
            self.assertEqual(int(row['assigned_kmers']), len(groups[taxid]))
            self.assertEqual(row['fasta_file'], f'taxon_fastas/taxid_{taxid}'+suffix)
            data = (output/row['fasta_file']).read_bytes()
            if compressed:
                self.assertTrue(data.startswith(b'\x1f\x8b'))
                data = gzip.decompress(data)  # Validates all members' CRCs/trailers.
            self.assertEqual(data, ''.join(groups[taxid]).encode())
        self.assertEqual({p.name for p in (output/'taxon_fastas').iterdir()},
                         {f'taxid_{t}'+suffix for t in groups})
        self.assertFalse((output/'taxon_fastas.partial').exists())
        self.assertFalse((output/'taxon_fastas.partial.tsv').exists())
        self.assertEqual(self.manifest(output)['write_separate_fastas'], 'yes')
        self.assertEqual(self.manifest(output)['separate_fasta_compression'], 'gzip' if compressed else 'none')
        self.assertEqual(int(self.manifest(output)['separate_fasta_files']), len(groups))
        self.assertTrue((output/'SUCCESS').exists())

    def test_separate_fastas_aliases_direct_assignments_and_optional_summary(self):
        baseline = self.root/'baseline'
        self.run_ok(self.blast_args(baseline, flag=None))
        self.assertFalse((baseline/'taxon_fastas').exists())
        for i, flag in enumerate(('-wsf', '--write-separate-fastas')):
            output = self.root/('split'+str(i))
            self.run_ok(self.blast_args(output, flag='-rd' if i else None)+[flag])
            self.check_separate(output, output/'markers.fasta')
            self.assertEqual((output/'markers.fasta').read_bytes(), (baseline/'markers.fasta').read_bytes())
            self.assertEqual((output/'markers.tsv').read_bytes(), (baseline/'markers.tsv').read_bytes())
            self.assertEqual((output/'descendant_summary.tsv').exists(), bool(i))
            if i:
                self.check_report(output)
            self.assertFalse((output/'taxon_fastas/taxid_9000003.fasta').exists())
        self.run_bad(self.blast_args(self.root/'duplicate')+['-wsf', '--write-separate-fastas'], 'Repeated option')

    def test_separate_fastas_from_portable_results_without_rd_or_references(self):
        baseline = self.root/'baseline'
        self.run_ok(self.blast_args(baseline, flag=None))
        before = {p: p.read_bytes() for p in baseline.iterdir() if p.is_file()}
        moved = self.root/'portable'; moved.mkdir()
        for name in ('markers.fasta', 'target_taxonomy.tsv', 'manifest.tsv', 'SUCCESS'):
            shutil.copyfile(baseline/name, moved/name)
        shutil.rmtree(self.db); self.fake.unlink()
        (self.root/'nodes.dmp').unlink(); (self.root/'names.dmp').unlink()
        output = self.root/'separate'
        result = subprocess.run([str(core.BINARY), '--write-separate-fastas', '-results', str(moved),
                                 '-out', str(output)], capture_output=True, text=True,
                                env=dict(os.environ, PATH=str(self.root/'empty-path')), timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('no reference scan', result.stderr)
        self.check_separate(output, moved/'markers.fasta')
        self.assertFalse((output/'descendant_summary.tsv').exists())
        for path, data in before.items():
            self.assertEqual(path.read_bytes(), data)
        combined = self.root/'combined'
        self.saved_report(moved, combined, extra=['-wsf'])
        self.check_separate(combined, moved/'markers.fasta')
        self.assertEqual((combined/'taxon_fastas.tsv').read_bytes(), (output/'taxon_fastas.tsv').read_bytes())

    def test_separate_fastas_empty_catalog(self):
        entries = [('a', 'ACGT'*30, '9000001', 'N/A'), ('b', 'ACGT'*30, '9606', 'N/A')]
        for compressed, flag in ((False, '-wsf'), (True, '-wsf-gz')):
            output = self.root/('empty'+str(compressed))
            self.run_ok(self.blast_args(output, entries)+[flag])
            self.check_separate(output, output/'markers.fasta', compressed)
            report = self.root/('empty-report'+str(compressed))
            self.saved_report(output, report, extra=[flag])
            self.check_separate(report, output/'markers.fasta', compressed)

    def test_separate_fastas_many_taxa_with_low_file_descriptor_limit(self):
        rng = random.Random(1200)
        nodes = [(1, 1), (10226, 1), (9606, 1)]+[(9100000+i, 10226) for i in range(128)]
        (self.root/'nodes.dmp').write_text(''.join(f'{t}\t|\t{p}\t|\n' for t, p in nodes))
        (self.root/'names.dmp').write_text(''.join(
            f'{t}\t|\tExample / same name; α\t|\t\t|\tscientific name\t|\n' for t, _ in nodes))
        shared = ''.join(rng.choices('ACGT', k=42))
        entries = [(str(t), 'N'.join([''.join(rng.choices('ACGT', k=42)) for _ in range(8)]+[shared]), str(t), 'N/A')
                   for t, _ in nodes[3:]]
        entries.append(('background', 'C'*42, '9606', 'N/A'))
        fasta = self.root/'many.fa'
        fasta.write_text(''.join('>'+e[0]+' taxid='+e[2]+'\n'+e[1]+'\n' for e in entries))
        def limit_files():
            _, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
            resource.setrlimit(resource.RLIMIT_NOFILE, (64 if hard == resource.RLIM_INFINITY else min(64, hard), hard))
        expected = {core.canonical(shared): 10226}
        for entry in entries[:-1]:
            for sequence in entry[1].split('N')[:-1]:
                expected[core.canonical(sequence)] = int(entry[2])
        for compressed, flag in ((False, '-wsf'), (True, '-wsf-gz')):
            output = self.root/('many'+str(compressed))
            args = self.common(output)+['-fasta', str(fasta), flag, '-rd']
            result = subprocess.run([str(core.BINARY)]+args, capture_output=True, text=True, preexec_fn=limit_files, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.check_separate(output, output/'markers.fasta', compressed)
            self.assertEqual({r['sequence']: int(r['specificity_taxid']) for r in self.table(output)}, expected)
            self.assertEqual(len(list((output/'taxon_fastas').iterdir())), 129)
            report = self.root/('many-report'+str(compressed))
            result = subprocess.run([str(core.BINARY), flag, '-results', str(output), '-out', str(report)],
                                    capture_output=True, text=True, preexec_fn=limit_files, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.check_separate(report, output/'markers.fasta', compressed)
            if compressed:
                max_members = 0
                paths = sorted((report/'taxon_fastas').iterdir())
                for path in paths:
                    rest = path.read_bytes(); members = 0
                    while rest:
                        decoder = zlib.decompressobj(31)
                        decoder.decompress(rest); decoder.flush()
                        self.assertTrue(decoder.eof, path.name)
                        rest = decoder.unused_data; members += 1
                    max_members = max(max_members, members)
                self.assertGreater(max_members, 1, 'Fixture must exercise gzip cache eviction/reopening')
                if core.GZIP_BINARY:
                    result = subprocess.run([core.GZIP_BINARY, '-t', '--']+[str(p) for p in paths], capture_output=True)
                    self.assertEqual(result.returncode, 0, result.stderr)

    def test_separate_fastas_checkpoint_persistence_and_output_rewrite(self):
        baseline = self.root/'baseline'
        self.run_ok(self.blast_args(baseline)+['-wsf'])
        for phase in ('target loading', 'background scan'):
            output = self.root/phase.replace(' ', '-')
            args = self.blast_args(output)+['-wsf']
            self.configure(expected_taxids=self.ids, **{
                'target_row_delay' if phase == 'target loading' else 'row_delay': 0.3})
            process, log = self.launch(args)
            self.wait_log(process, log, 'Saved '+phase+r' after [1-9][0-9]* records', regex=True)
            self.interrupt(process, log, signal.SIGKILL)
            self.configure(expected_taxids=self.ids)
            self.run_ok(['-resume', str(output)], binary=Path(os.environ.get('CLADEMERS_RESUME_BINARY', str(core.BINARY))))
            # The option must survive without being repeated, including on another architecture.
            self.compare(baseline, output)
            self.check_separate(output, output/'markers.fasta')
        plain = self.root/'plain'
        self.run_ok(self.blast_args(plain))
        (plain/'SUCCESS').unlink()
        for name in ('taxon_fastas', 'taxon_fastas.partial'):
            (plain/name).mkdir()
            (plain/name/'taxid_10226.fasta').write_text('stale partial data\n')
            (plain/name/'taxid_9999999.fasta').write_text('stale generated file\n')
        self.run_ok(['-resume', str(plain), '-wsf'])
        self.check_separate(plain, plain/'markers.fasta')
        self.assertEqual((plain/'taxon_fastas.tsv').read_bytes(), (baseline/'taxon_fastas.tsv').read_bytes())
        (plain/'SUCCESS').unlink()
        (plain/'taxon_fastas.partial').mkdir()
        sentinel = plain/'taxon_fastas.partial/notes.txt'; sentinel.write_text('preserve me')
        self.run_bad(['-resume', str(plain)], 'Unexpected file')
        self.assertEqual(sentinel.read_text(), 'preserve me')
        self.assertFalse((plain/'SUCCESS').exists())

    def test_separate_fastas_are_not_published_for_corrupt_saved_catalog(self):
        baseline = self.root/'baseline'
        self.run_ok(self.blast_args(baseline))
        path = baseline/'markers.fasta'
        path.write_bytes(path.read_bytes().replace(b'taxid=10227', b'taxid=10226', 1))
        for i, flag in enumerate(('-wsf', '-wsf-gz')):
            report = self.root/('damaged'+str(i))
            result = self.saved_report(baseline, report, success=False, extra=[flag])
            self.assertIn('fingerprint mismatch', result.stderr)
            self.assertFalse((report/'taxon_fastas').exists())
            self.assertFalse((report/'taxon_fastas.tsv').exists())
            self.assertFalse((report/'descendant_summary.tsv').exists())

    @unittest.skipUnless(os.environ.get('CLADEMERS_OLD_BINARY'), 'optional 0.11.0 binary')
    def test_separate_fastas_from_existing_011_results(self):
        baseline = self.root/'old-results'
        self.run_ok(self.blast_args(baseline, flag=None), binary=Path(os.environ['CLADEMERS_OLD_BINARY']))
        self.assertEqual(self.manifest(baseline)['clademers_version'], '0.11.0')
        for compressed, flag in ((False, '-wsf'), (True, '-wsf-gz')):
            report = self.root/('new-export'+str(compressed))
            self.saved_report(baseline, report, extra=[flag])
            self.check_separate(report, baseline/'markers.fasta', compressed)

    def test_gzip_separate_aliases_imply_wsf_and_preserve_catalogs(self):
        baseline = self.root/'plain'
        self.run_ok(self.blast_args(baseline)+['-wsf'])
        cases = (['-wsf-gz'], ['--write-separate-fastas-gz'], ['-wsf', '-wsf-gz'], ['-wsf-gz', '-wsf'])
        for i, flags in enumerate(cases):
            output = self.root/('gzip'+str(i))
            self.run_ok(self.blast_args(output, flag='-rd' if i else None)+flags)
            self.check_separate(output, output/'markers.fasta', True)
            self.assertEqual((output/'markers.fasta').read_bytes(), (baseline/'markers.fasta').read_bytes())
            self.assertEqual((output/'markers.tsv').read_bytes(), (baseline/'markers.tsv').read_bytes())
            self.assertEqual((output/'descendant_summary.tsv').exists(), bool(i))
            for path in (baseline/'taxon_fastas').iterdir():
                self.assertEqual(gzip.decompress((output/'taxon_fastas'/(path.name+'.gz')).read_bytes()), path.read_bytes())
        self.run_bad(self.blast_args(self.root/'duplicate-gzip')+['-wsf-gz', '--write-separate-fastas-gz'], 'Repeated option')

    def test_gzip_separate_portable_export_without_external_programs(self):
        baseline = self.root/'baseline'
        self.run_ok(self.blast_args(baseline, flag=None))
        moved = self.root/'portable'; moved.mkdir()
        for name in ('markers.fasta', 'target_taxonomy.tsv', 'manifest.tsv', 'SUCCESS'):
            shutil.copyfile(baseline/name, moved/name)
        before = {p: p.read_bytes() for p in moved.iterdir()}
        shutil.rmtree(self.db); self.fake.unlink()
        (self.root/'nodes.dmp').unlink(); (self.root/'names.dmp').unlink()
        output = self.root/'gz-export'
        result = subprocess.run([str(core.BINARY), '--write-separate-fastas-gz', '-results', str(moved),
                                 '-out', str(output)], capture_output=True, text=True,
                                env=dict(os.environ, PATH=str(self.root/'empty-path')), timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.check_separate(output, moved/'markers.fasta', True)
        self.assertFalse((output/'descendant_summary.tsv').exists())
        self.assertFalse((output/'commands.log').exists())
        self.assertFalse((output/'gzip.log').exists())
        for path, data in before.items():
            self.assertEqual(path.read_bytes(), data)

    def test_gzip_separate_checkpoint_persistence_and_format_switch(self):
        baseline = self.root/'baseline'
        self.run_ok(self.blast_args(baseline)+['-wsf-gz'])
        for phase in ('target loading', 'background scan'):
            output = self.root/('gzip-'+phase.replace(' ', '-'))
            args = self.blast_args(output)+['-wsf-gz']
            self.configure(expected_taxids=self.ids, **{
                'target_row_delay' if phase == 'target loading' else 'row_delay': 0.3})
            process, log = self.launch(args)
            self.wait_log(process, log, 'Saved '+phase+r' after [1-9][0-9]* records', regex=True)
            self.interrupt(process, log, signal.SIGKILL)
            self.configure(expected_taxids=self.ids)
            self.run_ok(['-resume', str(output)], binary=Path(os.environ.get('CLADEMERS_RESUME_BINARY', str(core.BINARY))))
            self.compare(baseline, output)
            self.check_separate(output, output/'markers.fasta', True)
        plain = self.root/'switch-to-gzip'
        self.run_ok(self.blast_args(plain)+['-wsf'])
        (plain/'SUCCESS').unlink()
        partial = plain/'taxon_fastas.partial'; partial.mkdir()
        (partial/'taxid_10226.fasta').write_text('stale plain data\n')
        (partial/'taxid_9000001.fasta.gz').write_bytes(b'stale truncated gzip')
        (plain/'taxon_fastas/taxid_9999999.fasta.gz').write_bytes(b'stale gzip')
        self.run_ok(['-resume', str(plain), '-wsf-gz'])
        self.check_separate(plain, plain/'markers.fasta', True)
        self.assertEqual((plain/'taxon_fastas.tsv').read_bytes(), (baseline/'taxon_fastas.tsv').read_bytes())

    def test_gzip_separate_close_failure_never_publishes_outputs(self):
        rng = random.Random(1300)
        sequence = 'N'.join(''.join(rng.choices('ACGT', k=42)) for _ in range(800))
        entries = [('target', sequence, '9000001', 'N/A'), ('background', 'C'*42, '9606', 'N/A')]
        baseline = self.root/'baseline'
        self.run_ok(self.blast_args(baseline, entries)+['-wsf-gz'])
        path = baseline/'taxon_fastas/taxid_9000001.fasta.gz'
        limit = path.stat().st_size - 1  # The final gzip trailer cannot finish.
        self.assertGreater(limit, 4096)  # Allows the small report manifest to be written.
        def limit_file_size():
            signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
            _, hard = resource.getrlimit(resource.RLIMIT_FSIZE)
            resource.setrlimit(resource.RLIMIT_FSIZE, (limit, hard))
        output = self.root/'failed-export'
        result = subprocess.run([str(core.BINARY), '-wsf-gz', '-results', str(baseline), '-out', str(output)],
                                capture_output=True, text=True, preexec_fn=limit_file_size, timeout=30)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Compressed FASTA close failed', result.stderr)
        self.assertFalse((output/'SUCCESS').exists())
        self.assertFalse((output/'taxon_fastas').exists())
        self.assertFalse((output/'taxon_fastas.tsv').exists())


if __name__ == '__main__':
    unittest.main(verbosity=2)
