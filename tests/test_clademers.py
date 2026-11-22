#!/usr/bin/env python3
"""Developer tests; Python is NOT a CladeMers runtime dependency.

The oracle works directly with DNA strings, independently of the C++ bit
encoding. Mocked BLAST exports test failure handling and difficult taxonomy.
Set BLAST_BIN to a real BLAST+ bin directory to also test a genuine version-5
database made locally by makeblastdb (no NT download or network required).
"""

import csv
import gzip
import json
import os
from pathlib import Path
import random
import shutil
import sqlite3
import struct
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
BINARY = Path(os.environ.get("CLADEMERS_BINARY", str(ROOT / "clademers"))).resolve()
SQLITE3_BINARY = os.environ.get("CLADEMERS_SQLITE3", shutil.which("sqlite3"))
GZIP_BINARY = os.environ.get("CLADEMERS_GZIP", shutil.which("gzip"))
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def canonical(sequence):
    return min(sequence, reverse_complement(sequence))


def kmers(sequence, k):
    sequence = sequence.upper()
    return {canonical(sequence[i:i+k]) for i in range(len(sequence)-k+1)
            if set(sequence[i:i+k]) <= set("ACGT")}


def row(accession, sequence, taxid="10226", leaves="N/A"):
    return f"{accession}\t{taxid}\t{leaves}\t{sequence}\n"


def common_provenance(records):
    # FASTA has record numbers but no separate BLAST leaf annotation.
    return [{k:v for k,v in r.items() if k not in ('source_record_1based','source_leaf_taxids')}
            for r in records]


def fnv1a64(data):
    value = 14695981039346656037
    for byte in data:
        value = ((value ^ byte) * 1099511628211) & ((1 << 64) - 1)
    return f'{value:016x}'


def write_taxonomy(directory):
    nodes = [(1, 1), (33208, 1), (10226, 33208), (10227, 10226),
             (9000001, 10227), (9606, 33208), (2, 1)]
    (directory / "nodes.dmp").write_text("".join(
        f"{node}\t|\t{parent}\t|\tno rank\t|\n" for node, parent in reversed(nodes)))
    (directory / "names.dmp").write_text("".join(
        f"{node}\t|\t{'Placozoa' if node == 10226 else 'taxon '+str(node)}\t|\t\t|\tscientific name\t|\n"
        for node, _ in nodes))


def write_blast_taxonomy(directory, dump_directory):
    """Independent fixture for NCBI's SQLite tree and big-endian taxdb format."""
    directory.mkdir(parents=True, exist_ok=True)
    nodes = [tuple(int(v.strip()) for v in line.split('|')[:2])
             for line in (dump_directory/'nodes.dmp').read_text().splitlines()]
    with sqlite3.connect(directory/'taxonomy4blast.sqlite3') as connection:
        connection.execute('CREATE TABLE TaxidInfo (taxid INTEGER PRIMARY KEY, parent INTEGER)')
        connection.executemany('INSERT INTO TaxidInfo VALUES (?,?)', nodes)
        connection.execute('CREATE INDEX TaxidInfoCompositeIdx_parent ON TaxidInfo(parent,taxid)')
    names = {}
    for line in (dump_directory/'names.dmp').read_text().splitlines():
        fields = [f.strip() for f in line.split('|')]
        if fields[3] == 'scientific name':
            names[int(fields[0])] = fields[1]
    index = bytearray(struct.pack('>6I', 0x8739, len(names), 0, 0, 0, 0))
    data = bytearray()
    for taxid, name in sorted(names.items()):
        index.extend(struct.pack('>2I', taxid, len(data)))
        data.extend((name + '\t\t\tEukaryota').encode())
    (directory/'taxdb.bti').write_bytes(index)
    (directory/'taxdb.btd').write_bytes(data)
    return directory


FAKE_BLAST = r'''
import json
from pathlib import Path
import sys
import time
args = sys.argv[1:]
if args == ['-version']:
    print('blastdbcmd: synthetic test fixture')
elif args == ['-help']:
    print('-no_taxid_expansion; %X means leaf-node taxids')
else:
    db_argument = args[args.index('-db')+1]
    assert db_argument.startswith('"') and db_argument.endswith('"')
    db = Path(db_argument[1:-1])
    config = json.loads((db/'config.json').read_text())
    if '-info' in args:
        count = db/'info_calls'
        n = int(count.read_text()) if count.exists() else 0
        count.write_text(str(n+1))
        if n and config.get('change_taxonomy_file'):
            with Path(config['change_taxonomy_file']).open('ab') as handle:
                handle.write(b'changed')
        print('test nucleotide database ' + ('changed' if config.get('changed') and n else 'stable'))
        if 'database_info' in config:
            print(config['database_info'])
    else:
        assert '-target_only' not in args
        target = '-taxidlist' in args
        assert args[args.index('-outfmt')+1] == ('' if target else '%o\t')+'%a\t%T\t%X\t%s'
        if target:
            assert '-no_taxid_expansion' in args
            ids = set(Path(args[args.index('-taxidlist')+1]).read_text().split())
            assert ids == set(config.get('expected_taxids', ['10226','10227','9000001']))
        rows = (db/('target.tsv' if target else 'all.tsv')).read_text().splitlines(keepends=True)
        for i,line in enumerate(rows):
            if target and i:
                time.sleep(config.get('target_row_delay',0))
            if not target:
                if i:
                    time.sleep(config.get('row_delay',0))
                line = str(config.get('oids',list(range(len(rows))))[i])+'\t'+line
            sys.stdout.write(line)
            sys.stdout.flush()
        sys.exit(config.get('target_exit' if target else 'background_exit', 0))
'''


class CladeMersTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="clademers-test-")
        self.root = Path(self.temp.name)
        write_taxonomy(self.root)
        self.fake = self.root / "fake blast'cmd.py"
        # sys.executable is absolute; the tested child is a real subprocess.
        self.fake.write_text("#!" + sys.executable + "\n" + FAKE_BLAST)
        self.fake.chmod(0o755)
        self.serial = 0

    def tearDown(self):
        self.temp.cleanup()

    def dump_arguments(self, nodes_compressed=False, names_compressed=False,
                       nodes_path=None, names_path=None):
        args, paths = [], {}
        for key, compressed, override in (('nodes',nodes_compressed,nodes_path),
                                          ('names',names_compressed,names_path)):
            path = Path(override) if override is not None else self.root/(key+'.dmp')
            if compressed and override is None:
                data = path.read_bytes()
                path = self.root/(f"{key} {self.serial} ' ; $(touch TAX_INJECTED).dmp.gz")
                path.write_bytes(gzip.compress(data,mtime=0))
            paths[key] = path
            args += ['-tax'+key,str(path)]
        return args, paths

    def check_dump_manifest(self, manifest, paths):
        self.assertEqual(manifest['taxonomy_inputs_match'],'yes')
        for key,path in paths.items():
            raw = path.read_bytes()
            compressed = raw.startswith(b'\x1f\x8b')
            data = gzip.decompress(raw) if compressed else raw
            self.assertEqual(manifest[key+'_compression'],'gzip' if compressed else 'none')
            self.assertEqual(manifest[key+'_fingerprint_scope'],'decompressed_content')
            self.assertEqual(manifest[key+'_fnv1a64'],fnv1a64(data))
            self.assertEqual(int(manifest[key+'_bytes']),len(data))

    def execute_fasta(self, content, mapping=None, k=42, threads=1, blasttax=None,
                      success=True, extra=None, compressed=False, suffix=None, gzip_reader=None,
                      mapping_compressed=False, map_suffix=None, nodes_compressed=False,
                      names_compressed=False, nodes_path=None, names_path=None):
        self.serial += 1
        suffix = suffix if suffix is not None else ('.fasta.gz' if compressed else '.fa')
        fasta = self.root / (f"references {self.serial} ' ; $(touch INJECTED)" + suffix)
        stored = content.encode() if isinstance(content, str) else content
        if compressed:
            stored = gzip.compress(stored, mtime=0)
        fasta.write_bytes(stored)
        is_gzip = stored.startswith(b'\x1f\x8b')
        output = self.root / f"fasta output {self.serial}"
        args = [str(BINARY), '-resume-mode', 'replay', '-fasta', str(fasta), '-taxid', '10226',
                '-k', str(k), '-t', str(threads), '-out', str(output)]
        reader = GZIP_BINARY if gzip_reader is None else gzip_reader
        dump_paths = {}
        if blasttax:
            args += ['-blasttax', str(blasttax), '-sqlite3', SQLITE3_BINARY]
        else:
            dump_args, dump_paths = self.dump_arguments(nodes_compressed,names_compressed,nodes_path,names_path)
            args += dump_args
        dump_gzips = sum(p.is_file() and p.read_bytes().startswith(b'\x1f\x8b') for p in dump_paths.values())
        map_is_gzip = False
        if mapping is not None:
            map_suffix = map_suffix if map_suffix is not None else ('.tsv.gz' if mapping_compressed else '.tsv')
            taxmap = self.root / (f"map {self.serial} ' ; $(touch INJECTED)" + map_suffix)
            map_stored = mapping.encode() if isinstance(mapping, str) else mapping
            if mapping_compressed:
                map_stored = gzip.compress(map_stored, mtime=0)
            map_is_gzip = map_stored.startswith(b'\x1f\x8b')
            taxmap.write_bytes(map_stored)
            args += ['-taxidmap', str(taxmap)]
        if (is_gzip or map_is_gzip or dump_gzips) and reader:
            args += ['-gzip', reader]
        # No blastdbcmd, shell, or implicit SQLite/gzip executable is available.
        env = dict(os.environ, PATH=str(self.root/'empty-path'))
        result = subprocess.run(args + (extra or []), capture_output=True, text=True,
                                cwd=self.root, env=env, timeout=60)
        self.last_result = result
        self.assertFalse((self.root/'INJECTED').exists())
        self.assertFalse((self.root/'TAX_INJECTED').exists())
        if not success:
            self.assertNotEqual(result.returncode, 0, result.stderr)
            self.assertFalse((output/'SUCCESS').exists())
            return result, output
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((output/'SUCCESS').is_file())
        with (output/'markers.tsv').open() as file:
            records = list(csv.DictReader(file, delimiter='\t'))
        with (output/'manifest.tsv').open() as file:
            manifest = {r['key']: r['value'] for r in csv.DictReader(file, delimiter='\t')}
        self.assertEqual(manifest['status'], 'complete')
        self.assertEqual(manifest['sequence_source'], 'fasta')
        self.assertEqual(manifest['fasta_taxid_source'], 'headers' if mapping is None else 'map')
        self.assertEqual(manifest['taxonomy_source'], 'blast' if blasttax else 'dmp')
        self.check_dump_manifest(manifest,dump_paths)
        self.assertEqual(manifest['fasta_passes_match'], 'yes')
        self.assertEqual(manifest['fasta_compression'], 'gzip' if is_gzip else 'none')
        self.assertEqual(manifest['fasta_fingerprint_scope'], 'decompressed_content')
        decoded = gzip.decompress(stored) if is_gzip else stored
        self.assertEqual(int(manifest['fasta_bytes']), len(decoded))
        self.assertEqual(int(manifest['retained_kmers']), len(records))
        # Check exact content-byte fingerprinting, including CRLF/EOF.
        self.assertEqual(manifest['fasta_fnv1a64'], fnv1a64(decoded))
        if mapping is not None:
            map_decoded = gzip.decompress(map_stored) if map_is_gzip else map_stored
            self.assertEqual(manifest['taxidmap_fnv1a64'], fnv1a64(map_decoded))
            self.assertEqual(int(manifest['taxidmap_bytes']), len(map_decoded))
            self.assertEqual(manifest['taxidmap_compression'], 'gzip' if map_is_gzip else 'none')
            self.assertEqual(manifest['taxidmap_fingerprint_scope'], 'decompressed_content')
            self.assertEqual(manifest['taxidmap_passes_match'], 'yes')
        self.assertEqual((output/'markers.fasta').read_text().splitlines()[1::2],
                         [r['sequence'] for r in records])
        self.assertFalse((output/'blastdbcmd.log').exists())
        self.assertFalse((output/'database_info.txt').exists())
        if is_gzip or map_is_gzip or dump_gzips:
            self.assertTrue((output/'gzip.log').exists())
            self.assertEqual(manifest['gzip'], reader)
            commands = (output/'commands.log').read_text().splitlines()
            self.assertEqual(sum('"-d" "-c" "--"' in line for line in commands),
                             2 * (int(is_gzip) + int(map_is_gzip) + dump_gzips))
        elif not blasttax:
            self.assertFalse((output/'commands.log').exists())
        return records, manifest, output

    @unittest.skipUnless(GZIP_BINARY, 'gzip CLI required for compressed dumps')
    def test_gzip_dump_all_text_combinations(self):
        content = (ROOT/'examples/references.fasta').read_text()
        mapping = (ROOT/'examples/taxid_map.tsv').read_text()
        reference, manifest, output = self.execute_fasta(content,mapping)
        ids = (output/'target_taxids.txt').read_bytes()
        for mask in range(16):
            with self.subTest(mask=mask):
                records, current, out = self.execute_fasta(content,mapping,
                    compressed=bool(mask&1),mapping_compressed=bool(mask&2),
                    nodes_compressed=bool(mask&4),names_compressed=bool(mask&8),threads=4)
                self.assertEqual(records,reference)
                self.assertEqual(current['target_name'],'Placozoa')
                self.assertEqual(current['nodes_fnv1a64'],manifest['nodes_fnv1a64'])
                self.assertEqual(current['names_fnv1a64'],manifest['names_fnv1a64'])
                self.assertEqual((out/'target_taxids.txt').read_bytes(),ids)
        headers = (ROOT/'examples/references_taxids.fasta').read_text()
        records, _, _ = self.execute_fasta(headers,compressed=True,
            nodes_compressed=True,names_compressed=True,threads=30)
        self.assertEqual(records,reference)

    @unittest.skipUnless(GZIP_BINARY, 'gzip CLI required for compressed dumps')
    def test_gzip_dump_blast_input_equivalence(self):
        target, background = row('target','ACGT'*20), row('bg','A'*80,'9606')
        reference, _, _ = self.execute(target,background)
        for nodes in (False,True):
            for names in (False,True):
                for workers in (1,30):
                    records, _, _ = self.execute(target,background,nodes_compressed=nodes,
                        names_compressed=names,threads=workers)
                    self.assertEqual(records,reference)
        # Explicit reader is permitted for BLAST+.dmp and unused when all inputs are plain.
        records, _, _ = self.execute(target,background,gzip_reader='/missing/gzip')
        self.assertEqual(records,reference)

    @unittest.skipUnless(GZIP_BINARY, 'gzip CLI required for compressed dumps')
    def test_gzip_dump_members_and_content_fingerprints(self):
        content = (ROOT/'examples/references_taxids.fasta').read_text()
        nodes = (self.root/'nodes.dmp').read_text().replace('\n','\r\n').rstrip().encode()
        names = ((self.root/'names.dmp').read_text()+
                 '2\t|\t'+'x'*70000+'\t|\t\t|\tsynonym\t|').replace('\n','\r\n').encode()
        (self.root/'nodes.dmp').write_bytes(nodes)
        (self.root/'names.dmp').write_bytes(names)
        reference, _, _ = self.execute_fasta(content)
        for key,data in (('nodes',nodes),('names',names)):
            for split in (1,len(data)//2,len(data)):
                path = self.root/(key+" ' ; $(touch TAX_INJECTED).payload")
                path.write_bytes(gzip.compress(data[:split],mtime=0)+gzip.compress(b'',mtime=0)+
                                 gzip.compress(data[split:],mtime=0))
                records, _, _ = self.execute_fasta(content,compressed=True,**{key+'_path':path})
                self.assertEqual(records,reference)

    @unittest.skipUnless(GZIP_BINARY, 'gzip CLI required for compressed dumps')
    def test_gzip_dump_invalid_inputs_and_taxonomy(self):
        content = (ROOT/'examples/references_taxids.fasta').read_text()
        for key in ('nodes','names'):
            raw = (self.root/(key+'.dmp')).read_bytes()
            packed = gzip.compress(raw,mtime=0)
            crc_bad = bytearray(packed)
            crc_bad[-8] ^= 1
            malformed = [b'',b'\x1f',b'\x1f\x8b',packed[:10],packed[:len(packed)//2],
                packed[:-1],packed[:-8],bytes(crc_bad),gzip.compress(b''),
                packed+gzip.compress(b'\n')[:-3],gzip.compress(b'bad fields\n'),raw]
            if key=='nodes':
                malformed += [gzip.compress(raw.replace(b'10227\t|\t10226',b'10227\t|\t9000001')),
                    gzip.compress(raw.replace(b'10226\t|\t33208',b'10226\t|\t444444')),
                    gzip.compress(raw+b'1\t|\t1\t|\tno rank\t|\n')]
            else:
                malformed += [gzip.compress(raw.replace(b'scientific name',b'synonym')),
                    gzip.compress(raw+b'10226\t|\tduplicate\t|\t\t|\tscientific name\t|\n')]
            path = self.root/(key+'.dmp.gz')
            for data in malformed:
                with self.subTest(key=key,bytes=len(data)):
                    path.write_bytes(data)
                    result,out=self.execute_fasta(content,success=False,**{key+'_path':path})
                    self.assertIn('ERROR:',result.stderr)
                    self.assertFalse((out/'markers.fasta').exists())

    @unittest.skipUnless(GZIP_BINARY, 'gzip CLI required for compressed dumps')
    def test_gzip_dump_reader_and_verification_failures(self):
        content = (ROOT/'examples/references_taxids.fasta').read_text()
        reader = self.root/"fake dump reader ' ; $(touch TAX_INJECTED).py"
        for key in ('nodes','names'):
            raw = (self.root/(key+'.dmp')).read_bytes()
            path = self.root/(key+'.payload')
            path.write_bytes(gzip.compress(raw,mtime=0))
            result,_=self.execute_fasta(content,gzip_reader='/missing/gzip',success=False,**{key+'_path':path})
            self.assertIn('Cannot start',result.stderr)
            for mode in ('load_exit','verify_exit','changed_output','changed_file','corrupted_file','changed_compression'):
                counter = self.root/(key+'_'+mode+'.count')
                path.write_bytes(gzip.compress(raw,mtime=0))
                reader.write_text('#!'+sys.executable+'\n'+
                    'import gzip,pathlib,sys\n'+
                    'assert sys.argv[1:4]==["-d","-c","--"]\n'+
                    'path=pathlib.Path(sys.argv[4])\nassert path.is_absolute()\n'+
                    'counter=pathlib.Path('+repr(str(counter))+')\n'+
                    'n=int(counter.read_text()) if counter.exists() else 0\n'+
                    'counter.write_text(str(n+1))\n'+
                    'data=gzip.decompress(path.read_bytes())\nmode='+repr(mode)+'\n'+
                    'if n and mode=="changed_output": data=data.replace(b"\\t",b" ",1)\n'+
                    'if not n and mode=="changed_file": path.write_bytes(gzip.compress(data.replace(b"\\t",b" ",1)))\n'+
                    'if not n and mode=="corrupted_file": path.write_bytes(path.read_bytes()[:-1])\n'+
                    'if not n and mode=="changed_compression": path.write_bytes(data)\n'+
                    'sys.stdout.buffer.write(data)\n'+
                    'sys.exit(9 if (mode=="load_exit" or (n and mode=="verify_exit")) else 0)\n')
                reader.chmod(0o755)
                result,out=self.execute_fasta(content,gzip_reader=str(reader),success=False,**{key+'_path':path})
                self.assertIn('Taxonomy file changed' if mode.startswith('changed') else 'did not complete',result.stderr)
                self.assertIn('status\tfailed',(out/'manifest.tsv').read_text())
                self.assertFalse((out/'markers.fasta').exists())

    @unittest.skipUnless(GZIP_BINARY, 'gzip CLI required for compressed maps')
    def test_gzip_map_all_input_combinations(self):
        content = (ROOT/'examples/references.fasta').read_text()
        mapping = (ROOT/'examples/taxid_map.tsv').read_text()
        taxonomy = write_blast_taxonomy(self.root/'blast taxonomy', self.root) if SQLITE3_BINARY else None
        reference, manifest, _ = self.execute_fasta(content, mapping)
        for tax in ([None,taxonomy] if taxonomy else [None]):
            for fasta_gzip in (False,True):
                for map_gzip in (False,True):
                    for workers in (1,30):
                        records, current, _ = self.execute_fasta(content, mapping, blasttax=tax,
                            compressed=fasta_gzip, mapping_compressed=map_gzip, threads=workers)
                        self.assertEqual(records,reference)
                        self.assertEqual(current['taxidmap_fnv1a64'],manifest['taxidmap_fnv1a64'])

    @unittest.skipUnless(GZIP_BINARY, 'gzip CLI required for compressed maps')
    def test_gzip_map_members_wrapping_and_exact_content(self):
        content = (ROOT/'examples/references.fasta').read_text()
        # Long ignored line, comments, CRLF, optional header, no final newline.
        mapping = ('#'+'x'*70000+'\r\nsequence_id\ttaxids\r\n\r\n'+
                   (ROOT/'examples/taxid_map.tsv').read_text().replace('\n','\r\n').rstrip()).encode()
        reference, manifest, _ = self.execute_fasta(content, mapping)
        for split in (1, 65536, len(mapping)-2, len(mapping)):
            packed = gzip.compress(mapping[:split],mtime=0)+gzip.compress(b'',mtime=0)
            packed += gzip.compress(mapping[split:],mtime=0)
            records, current, _ = self.execute_fasta(content, packed, map_suffix='.data', compressed=True)
            self.assertEqual(records,reference)
            self.assertEqual(current['taxidmap_fnv1a64'],manifest['taxidmap_fnv1a64'])

    @unittest.skipUnless(GZIP_BINARY, 'gzip CLI required for compressed maps')
    def test_gzip_map_union_missing_and_conflicting_taxids(self):
        rng = random.Random(605)
        seqs = [''.join(rng.choices('ACGT',k=42)) for _ in range(6)]
        content = ''.join(f'>t{i}\n{s}\n' for i,s in enumerate(seqs))
        content += ''.join(f'>{name}\n{s}\n' for name,s in zip(
            ['mixed','unknown','missing taxid=10227','version.1','pipe|id'],seqs))
        mapping = ''.join(f't{i} 10227\n' for i in range(6))
        mapping += 'mixed 10227\nmixed 9606\nunknown N/A\nversion 10227\npipe 10227\n'
        records, _, _ = self.execute_fasta(content, mapping, compressed=True, mapping_compressed=True)
        self.assertEqual([r['sequence'] for r in records],[canonical(seqs[-1])])
        agrees = '>t taxids=9000001,10227\n'+'ACGT'*20+'\n>b taxid=2\n'+'N'*80
        self.execute_fasta(agrees, 't 10227\nt 9000001\nb 2', mapping_compressed=True)
        result, _ = self.execute_fasta(agrees, 't 10227\nb 2', mapping_compressed=True,success=False)
        self.assertIn('Header and mapping taxids disagree',result.stderr)

    @unittest.skipUnless(GZIP_BINARY, 'gzip CLI required for compressed maps')
    def test_gzip_map_corruption_and_format_errors(self):
        content = (ROOT/'examples/references.fasta').read_text()
        mapping = (ROOT/'examples/taxid_map.tsv').read_bytes()
        packed = gzip.compress(mapping,mtime=0)
        crc_bad = bytearray(packed)
        crc_bad[-8] ^= 1
        malformed = [b'',b'\x1f',b'\x1f\x8b',packed[:10],packed[:len(packed)//2],packed[:-1],
                     packed[:-8],bytes(crc_bad),gzip.compress(b''),packed+gzip.compress(b'\n')[:-3],
                     gzip.compress(b'sequence_id taxid\n'),gzip.compress(b't 10227 extra\n'),
                     gzip.compress(b't 10227,,9606\n'),gzip.compress(b't 4294967296\n')]
        for data in malformed:
            with self.subTest(bytes=len(data)):
                result, out = self.execute_fasta(content, data, map_suffix='.tsv.gz',success=False)
                self.assertIn('ERROR:',result.stderr)
                self.assertFalse((out/'markers.fasta').exists())
        result, _ = self.execute_fasta(content, mapping, map_suffix='.tsv.gz',success=False)
        self.assertIn('no gzip header',result.stderr)

    @unittest.skipUnless(GZIP_BINARY, 'gzip CLI required for compressed maps')
    def test_gzip_map_reader_and_verification_failures(self):
        content = (ROOT/'examples/references.fasta').read_text()
        mapping = (ROOT/'examples/taxid_map.tsv').read_text()
        result, _ = self.execute_fasta(content,mapping,mapping_compressed=True,
                                      gzip_reader='/missing/gzip',success=False)
        self.assertIn('Cannot start',result.stderr)
        reader = self.root/"fake map reader ' ; $(touch INJECTED).py"
        for mode in ('load_exit','verify_exit','changed_taxid','changed_newline','changed_file','corrupted_file'):
            counter = self.root/(mode+'.count')
            reader.write_text('#!'+sys.executable+'\n'+
                'import gzip,pathlib,sys\n'+
                'assert sys.argv[1:4] == ["-d","-c","--"]\n'+
                'path=pathlib.Path(sys.argv[4])\nassert path.is_absolute()\n'+
                'counter=pathlib.Path('+repr(str(counter))+')\n'+
                'n=int(counter.read_text()) if counter.exists() else 0\n'+
                'counter.write_text(str(n+1))\n'+
                'data=gzip.decompress(path.read_bytes())\nmode='+repr(mode)+'\n'+
                'if n and mode=="changed_taxid": data=data.replace(b"10227",b"10226")\n'+
                'if n and mode=="changed_newline": data=data.rstrip(b"\\n")\n'+
                'if not n and mode=="changed_file": path.write_bytes(gzip.compress(data.replace(b"10227",b"10226")))\n'+
                'if not n and mode=="corrupted_file": path.write_bytes(path.read_bytes()[:-1])\n'+
                'sys.stdout.buffer.write(data)\n'+
                'sys.exit(9 if (mode=="load_exit" or (n and mode=="verify_exit")) else 0)\n')
            reader.chmod(0o755)
            result,out=self.execute_fasta(content,mapping,mapping_compressed=True,
                                          gzip_reader=str(reader),success=False)
            self.assertIn('mapping file changed' if mode.startswith('changed') else 'did not complete',result.stderr)
            self.assertIn('status\tfailed',(out/'manifest.tsv').read_text())
            self.assertFalse((out/'markers.fasta').exists())
        # With both inputs plain, even an unavailable gzip reader is unused.
        self.execute_fasta(content,mapping,extra=['-gzip','/missing/gzip'])

    @unittest.skipUnless(GZIP_BINARY, 'gzip CLI required for compressed FASTA')
    def test_gzip_header_map_taxonomy_and_worker_equivalence(self):
        taxonomy = write_blast_taxonomy(self.root/'blast taxonomy', self.root) if SQLITE3_BINARY else None
        for tax in ([None,taxonomy] if taxonomy else [None]):
            for mapped in (False,True):
                content = (ROOT/'examples'/('references.fasta' if mapped else 'references_taxids.fasta')).read_text()
                mapping = (ROOT/'examples/taxid_map.tsv').read_text() if mapped else None
                plain, manifest, _ = self.execute_fasta(content, mapping, blasttax=tax)
                for workers in (1,30):
                    records, compressed_manifest, _ = self.execute_fasta(content, mapping, blasttax=tax,
                        compressed=True, threads=workers)
                    self.assertEqual(records, plain)
                    self.assertEqual(manifest['fasta_fnv1a64'], compressed_manifest['fasta_fnv1a64'])
                    self.assertEqual(compressed_manifest['retained_kmers'], '1')

    @unittest.skipUnless(GZIP_BINARY, 'gzip CLI required for compressed FASTA')
    def test_gzip_string_oracle_and_chunk_boundaries(self):
        rng = random.Random(746)
        sequence = ''.join(rng.choices('ACGT',k=200))
        target = 'NN' + sequence[:100].lower() + 'r-U' + sequence[100:]
        # A long line crosses the process buffer and scanner's 1-MiB boundary.
        background = 'N'*(1024*1024-20) + reverse_complement(sequence[40:160]) + 'N'*80
        content = f'>t|kraken:taxid|10227\r\n{target[:60]} \t\r\n{target[60:]}\r\n>b\r\n{background}'
        for k in (1,2,7,31,32,33,42,63,64):
            with self.subTest(k=k):
                records, _, _ = self.execute_fasta(content, k=k, compressed=True, threads=4)
                self.assertEqual({r['sequence'] for r in records}, kmers(target,k)-kmers(background,k))
                for r in records:
                    start = int(r['start_1based'])-1
                    window = target[start:start+k].upper()
                    self.assertEqual(r['sequence'], reverse_complement(window) if r['strand']=='-' else window)

    @unittest.skipUnless(GZIP_BINARY, 'gzip CLI required for compressed FASTA')
    def test_gzip_magic_concatenated_members_and_final_line(self):
        content = (ROOT/'examples/references_taxids.fasta').read_bytes().rstrip(b'\n')
        plain, manifest, _ = self.execute_fasta(content)
        # A member boundary can occur inside a header, inside DNA, or at EOF.
        for split in (1, 10, len(content)//2, len(content)):
            packed = gzip.compress(content[:split],mtime=0) + gzip.compress(b'',mtime=0)
            packed += gzip.compress(content[split:],mtime=0)
            records, compressed_manifest, _ = self.execute_fasta(packed, suffix='.data')
            self.assertEqual(records,plain)
            self.assertEqual(compressed_manifest['fasta_fnv1a64'], manifest['fasta_fnv1a64'])
        records, _, _ = self.execute_fasta(content, compressed=True, suffix='.fa')
        self.assertEqual(records,plain) # magic detection, not a filename requirement

    @unittest.skipUnless(GZIP_BINARY, 'gzip CLI required for compressed FASTA')
    def test_gzip_corrupt_truncated_and_empty_inputs(self):
        content = (ROOT/'examples/references_taxids.fasta').read_bytes()
        packed = gzip.compress(content,mtime=0)
        crc_bad = bytearray(packed)
        crc_bad[-8] ^= 1
        malformed = [b'', b'\x1f', b'\x1f\x8b', packed[:10], packed[:len(packed)//2],
                     packed[:-1], packed[:-8], bytes(crc_bad), gzip.compress(b''),
                     packed + gzip.compress(b'\n')[:-3], gzip.compress(b'not fasta\n')]
        for data in malformed:
            with self.subTest(bytes=len(data)):
                result, out = self.execute_fasta(data, suffix='.fasta.gz', success=False)
                self.assertIn('ERROR:', result.stderr)
                self.assertFalse((out/'markers.fasta').exists())
        # A misleading .gz suffix is an error; no silent plain-text fallback.
        result, _ = self.execute_fasta(content, suffix='.fasta.gz', success=False)
        self.assertIn('no gzip header', result.stderr)

    @unittest.skipUnless(GZIP_BINARY, 'gzip CLI required for compressed FASTA')
    def test_gzip_reader_failures_and_changed_second_pass(self):
        content = (ROOT/'examples/references_taxids.fasta').read_text()
        result, _ = self.execute_fasta(content, compressed=True, gzip_reader='/missing/gzip', success=False)
        self.assertIn('Cannot start', result.stderr)
        # A controlled producer verifies that EOF/output alone is never success.
        reader = self.root/"fake gzip ' ; $(touch INJECTED).py"
        for mode in ('first_exit', 'second_exit', 'changed_content', 'changed_newline'):
            counter = self.root/(mode+'.count')
            reader.write_text('#!'+sys.executable+'\n'+
                'import pathlib,sys\n'+
                'assert sys.argv[1:4] == ["-d","-c","--"]\n'+
                'assert pathlib.Path(sys.argv[4]).is_absolute()\n'+
                'counter=pathlib.Path('+repr(str(counter))+')\n'+
                'n=int(counter.read_text()) if counter.exists() else 0\n'+
                'counter.write_text(str(n+1))\n'+
                'content='+repr(content)+'\n'+
                'mode='+repr(mode)+'\n'+
                'if n and mode=="changed_content": content=content.replace("synthetic","Synthetic",1)\n'+
                'if n and mode=="changed_newline": content=content.rstrip("\\n")\n'+
                'sys.stdout.write(content)\n'+
                'sys.exit(9 if (mode=="first_exit" or (n and mode=="second_exit")) else 0)\n')
            reader.chmod(0o755)
            result, out = self.execute_fasta(content, compressed=True, gzip_reader=str(reader), success=False)
            self.assertIn('changed between passes' if mode.startswith('changed') else 'did not complete', result.stderr)
            self.assertIn('status\tfailed', (out/'manifest.tsv').read_text())
            self.assertFalse((out/'markers.fasta').exists())
        # Plain FASTA does not try to execute a gzip reader.
        self.execute_fasta(content, extra=['-gzip','/missing/gzip'])
        self.execute('', '', extra=['-gzip','/missing/gzip'], success=False)

    def test_fasta_headers_oracle_boundaries_wrapping_and_provenance(self):
        rng = random.Random(3041)
        seqs = [''.join(rng.choices('ACGT', k=180)) for _ in range(5)]
        seqs[0] = 'n-' + seqs[0][:80].lower() + 'Uryswkmbdhv' + seqs[0][80:]
        headers = ['same_id taxid=10226', 'same_id [taxid=10227]',
                   'three taxids=10227;9000001 taxids=9000001,10227',
                   'four|kraken:taxid|9000001 desc', 'kraken:taxid|10227|five [taxid=10227]']
        backgrounds = [seqs[1][20:100], reverse_complement(seqs[2][40:120]), 'N'*80]
        # Blank lines, tabs/spaces in sequence, CRLF, no final newline.
        content = '\r\n'.join('>'+header+'\r\n'+seq[:50]+' \t\r\n'+seq[50:]
                                for header, seq in zip(headers, seqs))
        content += ''.join(f'\r\n>bg{i} taxid=9606\r\n{s}' for i,s in enumerate(backgrounds))
        for k in (1, 2, 7, 31, 32, 33, 42, 63, 64):
            with self.subTest(k=k):
                expected = set().union(*(kmers(s,k) for s in seqs))
                expected -= set().union(*(kmers(s,k) for s in backgrounds))
                serial, manifest, _ = self.execute_fasta(content, k=k)
                parallel, _, _ = self.execute_fasta(content, k=k, threads=4)
                self.assertEqual(serial, parallel)
                self.assertEqual({r['sequence'] for r in serial}, expected)
                self.assertEqual(manifest['fasta_records'], '8')
                for r in serial:
                    index = int(r['source_record_1based'])-1
                    start = int(r['start_1based'])-1
                    window = seqs[index][start:start+k].upper()
                    self.assertEqual(r['sequence'], reverse_complement(window) if r['strand']=='-' else window)
                    self.assertEqual(r['source_accession'], headers[index].split()[0])
                    self.assertEqual(r['source_leaf_taxids'], 'N/A')

    def test_fasta_mixed_unknown_and_unrecognized_headers_are_background(self):
        rng = random.Random(778)
        seqs = [''.join(rng.choices('ACGT', k=42)) for _ in range(11)]
        content = ''.join(f'>t{i} taxid=10227\n{s}\n' for i,s in enumerate(seqs))
        headers = ['b0 taxids=10227,9606', 'b1 taxids=10227;0', 'b2 taxid=N/A',
                   'b3 taxid=0', 'b4 taxid=33208', 'b5 taxid=4000000000',
                   'b6 Placozoa 10227', '10227', 'b8 taxids=10227,N/A',
                   'b9 [taxids=10227,9000001]']
        content += ''.join(f'>{h}\n{s}\n' for h,s in zip(headers,seqs))
        records, manifest, _ = self.execute_fasta(content, threads=30)
        self.assertEqual({r['sequence'] for r in records}, {canonical(s) for s in seqs[9:]})
        self.assertEqual(manifest['fasta_unassigned_taxid_records'], '6')

    def test_fasta_mapping_union_exact_ids_and_authority(self):
        rng = random.Random(889)
        seqs = [''.join(rng.choices('ACGT',k=42)) for _ in range(7)]
        content = ''.join(f'>t{i}\n{s}\n' for i,s in enumerate(seqs))
        content += ''.join(f'>{h}\n{s}\n' for h,s in zip(
            ['mixed', 'unknown', 'missing taxid=10227', 'version.1', 'pipe|id', 'alltarget'], seqs))
        mapping = '# two-column map\r\nsequence_id\ttaxids\r\n\r\n'
        mapping += ''.join(f't{i}\t10227\r\n' for i in range(len(seqs)))
        mapping += ('mixed\t10227\r\nmixed\t9606\r\nunknown\tN/A\r\n'
                    'version\t10227\r\npipe\t10227\r\n'
                    'alltarget\t10227,9000001\r\nalltarget\t10227\r\nunused\t2')
        records, manifest, _ = self.execute_fasta(content, mapping, threads=4)
        self.assertEqual({r['sequence'] for r in records}, {canonical(s) for s in seqs[5:]})
        self.assertEqual(manifest['fasta_unassigned_taxid_records'], '4')
        self.assertEqual(manifest['taxidmap_ids'], '13')
        self.assertIn('taxidmap_fnv1a64', manifest)
        # Repeated records may share an ID; provenance remains unambiguous.
        records, _, _ = self.execute_fasta('>same\n'+seqs[0]+'\n>same\n'+seqs[1]+
            '\n>bg\n'+'N'*42, 'same 10227\nbg 2\n')
        self.assertEqual({r['source_record_1based'] for r in records}, {'1','2'})

    def test_fasta_header_and_mapping_agreement(self):
        content = '>t taxids=9000001,10227\n'+'ACGT'*20+'\n>b taxid=2\n'+'N'*80
        records, _, _ = self.execute_fasta(content, 't 10227\nt 9000001\nb 2\n')
        self.assertTrue(records)
        for mapping in ('t 10227\nb 2', 't 10227,9606\nb 2', 't 0\nb 2'):
            result, _ = self.execute_fasta(content, mapping, success=False)
            self.assertIn('Header and mapping taxids disagree', result.stderr)

    def test_fasta_rejects_malformed_or_uninformative_input(self):
        valid = '>t taxid=10227\n'+'ACGT'*20+'\n>b taxid=2\n'+'N'*80
        malformed = ['', 'ACGT\n', '>\nACGT\n', '>t taxid=10227\n>b\nAAAA',
                     valid+'\n>empty', valid+'\nX', valid+'\n1', valid+'\n*',
                     '>t taxid=10227\x01b taxid=2\nACGT',
                     '>t taxid=10227\nACGT'*1, '>t taxid=10227\n'+'A'*80,
                     '>t\n'+'A'*80+'\n>b\n'+'C'*80,
                     b'\x1f\x8bcompressed-data']
        for content in malformed:
            with self.subTest(content=repr(content)[:50]):
                self.execute_fasta(content, success=False)
        for annotation in ('taxid=', 'taxid=-1', 'taxid=4294967296', 'taxid=abc',
                           'taxids=10227,', 'taxids=10227,,9606', 'taxids=10227; 9606',
                           '[taxid=10227', 'taxid=10227 taxid=9606',
                           'taxid=10227 [taxids=10227,9606]'):
            with self.subTest(annotation=annotation):
                result, _ = self.execute_fasta(valid.replace('taxid=10227', annotation), success=False)
                self.assertIn('ERROR:', result.stderr)
        for identifier in ('t|kraken:taxid|', 'kraken:taxid|invalid', 't|kraken:taxid|10227,9606',
                           't|kraken:taxid|10227 taxid=9606'):
            self.execute_fasta(valid.replace('t taxid=10227', identifier), success=False)
        for mapping in ('', '# comment\nsequence_id taxid\n', 't', 't 10227 extra',
                        't -1', 't 4294967296', 't 10227,', 't 10227,,2', 't bad',
                        't\x00 10227', b'\x1f\x8bcompressed-data'):
            with self.subTest(mapping=repr(mapping)):
                self.execute_fasta(valid, mapping=mapping, success=False)

    def test_fasta_cli_exclusivity_and_missing_files(self):
        valid = '>t taxid=10227\n'+'ACGT'*20+'\n>b taxid=2\n'+'N'*80
        for extra in (['-blastdb','db'], ['-blastdbcmd','anything'], ['-fasta','missing'],
                      ['-taxidmap','missing'], ['-sqlite3','anything']):
            self.execute_fasta(valid, extra=extra, success=False)
        self.execute('', '', extra=['-taxidmap','missing'], success=False)
        for fasta in ('-', str(self.root), str(self.root/'missing')):
            output = self.root/'missing-result'
            result = subprocess.run([str(BINARY),'-fasta',fasta,'-taxnodes',str(self.root/'nodes.dmp'),
                '-taxnames',str(self.root/'names.dmp'),'-taxid','10226','-out',str(output)],
                capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertFalse((output/'SUCCESS').exists())
            if output.exists():
                shutil.rmtree(output)

    @unittest.skipUnless(SQLITE3_BINARY, 'sqlite3 CLI required for BLAST taxonomy')
    def test_fasta_both_taxonomies_and_annotation_sources(self):
        taxonomy = write_blast_taxonomy(self.root/'blast taxonomy', self.root)
        plain = (ROOT/'examples/references.fasta').read_text()
        mapping = (ROOT/'examples/taxid_map.tsv').read_text()
        lookup = dict(line.split() for line in mapping.splitlines())
        headers = '\n'.join(line+' taxid='+lookup[line[1:].split()[0]] if line.startswith('>') else line
                            for line in plain.splitlines())+'\n'
        reference = None
        for blasttax in (None, taxonomy):
            for mapped in (False,True):
                for workers in (1,30):
                    records, manifest, _ = self.execute_fasta(plain if mapped else headers,
                        mapping if mapped else None, threads=workers, blasttax=blasttax)
                    self.assertEqual(manifest['candidate_kmers'],'2')
                    self.assertEqual(manifest['retained_kmers'],'1')
                    if reference is None:
                        reference = records
                    self.assertEqual(records,reference)

    def execute(self, target, background, k=42, threads=1, config=None,
                success=True, extra=None, database=None, blast=None, blasttax=None,
                target_taxid=10226, sqlite3_reader=None, nodes_compressed=False,
                names_compressed=False, nodes_path=None, names_path=None, gzip_reader=None):
        self.serial += 1
        if database is None:
            database = self.root / (f"db {self.serial} ' ; $(touch INJECTED)")
            database.mkdir()
            (database / "target.tsv").write_text(target)
            (database / "all.tsv").write_text(background)
            (database / "config.json").write_text(json.dumps(config or {}))
        output = self.root / f"output {self.serial}"
        dump_paths = {}
        if blasttax:
            taxonomy_args = ['-blasttax',str(blasttax)]
        else:
            taxonomy_args, dump_paths = self.dump_arguments(nodes_compressed,names_compressed,nodes_path,names_path)
        dump_gzips = sum(p.is_file() and p.read_bytes().startswith(b'\x1f\x8b') for p in dump_paths.values())
        if gzip_reader is not None or (dump_gzips and GZIP_BINARY):
            taxonomy_args += ['-gzip', gzip_reader if gzip_reader is not None else GZIP_BINARY]
        reader = sqlite3_reader or SQLITE3_BINARY
        if blasttax and reader:
            taxonomy_args += ["-sqlite3", reader]
        command = [str(BINARY), '-resume-mode', 'replay'] + taxonomy_args + ["-taxid", str(target_taxid),
                   "-k", str(k), "-t", str(threads), "-blastdb", str(database),
                   "-out", str(output), "-blastdbcmd", str(blast or self.fake)]
        if extra:
            command += extra
        result = subprocess.run(command, text=True, capture_output=True, cwd=self.root, timeout=60)
        self.last_result = result
        self.assertFalse((self.root/"INJECTED").exists())
        self.assertFalse((self.root/"TAX_INJECTED").exists())
        if success:
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue((output/"SUCCESS").is_file())
            with (output/"markers.tsv").open() as file:
                records = list(csv.DictReader(file, delimiter="\t"))
            with (output/"manifest.tsv").open() as file:
                manifest = {r["key"]: r["value"] for r in csv.DictReader(file, delimiter="\t")}
            self.assertEqual(manifest["status"], "complete")
            self.assertEqual(int(manifest["retained_kmers"]), len(records))
            ids = (config or {}).get('expected_taxids', ['10226', '10227', '9000001'])
            self.assertEqual((output/"target_taxids.txt").read_text(),
                             ''.join(str(t)+'\n' for t in sorted(map(int, ids))))
            self.assertEqual(manifest['taxonomy_source'], 'blast' if blasttax else 'dmp')
            self.check_dump_manifest(manifest,dump_paths)
            fasta = (output/"markers.fasta").read_text().splitlines()
            self.assertEqual(fasta[1::2], [r["sequence"] for r in records])
            self.assertEqual(fasta[::2], [">"+r["marker_id"]+" taxid="+r["specificity_taxid"] for r in records])
            self.assertEqual([r["sequence"] for r in records], sorted(r["sequence"] for r in records))
            return records, manifest, output
        self.assertNotEqual(result.returncode, 0, result.stderr)
        self.assertFalse((output/"SUCCESS").exists())
        return result, output

    def test_progress_plain_fasta_and_output(self):
        content = (ROOT/'examples/references.fasta').read_text()
        mapping = (ROOT/'examples/taxid_map.tsv').read_text()
        records, manifest, _ = self.execute_fasta(content, mapping, extra=['-progress-seconds','1'])
        log = self.last_result.stderr
        self.assertEqual(self.last_result.stdout, '')
        self.assertEqual(manifest['progress_seconds'], '1')
        self.assertIn('named operation', manifest['progress_scope'])
        for label in ('FASTA target pass', 'FASTA background input'):
            self.assertIn(label+': 0.0%', log)
            self.assertRegex(log, label+r': 100\.0% \| elapsed \S+ \| ETA 00:00:00 \| '+
                             f'{len(content)}/{len(content)} content bytes')
        self.assertIn('Read taxid map: 0.0%', log)
        self.assertIn('Collect retained markers: 100.0%', log)
        self.assertIn(f'{len(records)}/{len(records)} markers', log)
        self.assertIn('Write retained markers: 100.0%', log)
        self.assertNotIn('\r', log)

    @unittest.skipUnless(GZIP_BINARY, 'gzip CLI required for progress test')
    def test_progress_gzip_unknown_first_pass_known_second_pass(self):
        content = (ROOT/'examples/references.fasta').read_bytes()
        mapping = (ROOT/'examples/taxid_map.tsv').read_bytes()
        # Concatenated members have no single footer containing the full size.
        packed = gzip.compress(content[:43])+gzip.compress(content[43:])
        self.execute_fasta(packed, mapping, mapping_compressed=True,
                           nodes_compressed=True, names_compressed=True)
        log = self.last_result.stderr
        self.assertRegex(log, r'FASTA target pass: percent unavailable.*ETA unavailable \(total unknown\)')
        self.assertRegex(log, r'FASTA background input: 0\.0%.*0/'+str(len(content))+r' content bytes')
        self.assertRegex(log, r'Read taxid map: percent unavailable.*ETA unavailable')
        for key in ('nodes','names','map'):
            self.assertRegex(log, r'Recheck '+key+r'[^\n]*: 0\.0%.*0/\d+ content bytes')
        self.assertEqual(self.last_result.stdout, '')

    def test_progress_blast_eta_and_duplicate_oids(self):
        target = row('t', 'ACGT'*100)
        background = target+row('same_oid', 'ACGT'*100, '9606')+row('other', 'A'*400, '9606')
        records, manifest, _ = self.execute(target, background, extra=['-progress-seconds','1'],
            config={'database_info':'\t2 sequences; 800 total bases', 'oids':[0,0,5], 'row_delay':1.1})
        self.assertEqual(records, [])  # Both deflines still participate in subtraction.
        self.assertEqual(manifest['blast_progress_total_bases'], '800')
        self.assertEqual(manifest['all_defline_rows'], '3')
        log = self.last_result.stderr
        self.assertRegex(log, r'BLAST background scan: 50\.0%.*ETA ~\d\d:\d\d:\d\d.*400/800 sequence bases')
        self.assertRegex(log, r'BLAST background scan: 100\.0%.*ETA 00:00:00.*800/800 sequence bases')
        self.assertIn('BLAST target loading: percent unavailable', log)
        self.assertIn('Write retained markers: 100.0%', log)
        self.assertIn('0/0 markers', log)
        lines = [line for line in log.splitlines() if 'BLAST background scan:' in line]
        self.assertEqual(sum(': 100.0%' in line for line in lines), 1)
        self.assertNotIn('nan', '\n'.join(lines).lower())
        self.assertNotIn('inf', '\n'.join(lines).lower())

    def test_progress_blast_metadata_fallbacks(self):
        target, background = row('t','ACGT'*20), row('b','A'*80,'9606')
        for hint in ('', 'unrecognized metadata', '1 sequences; 0 total bases',
                     '1 sequences; 999999999999999999999999 total bases'):
            with self.subTest(hint=hint):
                self.execute(target, background, config={'database_info':hint})
                self.assertIn('BLAST background scan: percent unavailable', self.last_result.stderr)
        _, manifest, _ = self.execute(target, background,
            config={'database_info':'1,234 sequences; 12,345,678 total bases'})
        self.assertEqual(manifest['blast_progress_total_bases'], '12345678')

    def test_progress_unreliable_totals_and_failed_export(self):
        target = row('t', 'ACGT'*20)
        background = row('a','A'*80,'9606')+row('b','C'*80,'9606')
        for hint,oids in (('2 sequences; 100 total bases',[0,1]),
                          ('2 sequences; 160 total bases',[5,3])):
            self.execute(target, background, extra=['-progress-seconds','1'],
                config={'database_info':hint, 'oids':oids, 'row_delay':1.1})
            log = self.last_result.stderr
            self.assertRegex(log, r'BLAST background scan: percent unavailable.*ETA unavailable')
        result, _ = self.execute(target, background, success=False,
            config={'database_info':'2 sequences; 160 total bases', 'background_exit':9})
        self.assertNotIn('BLAST background scan: 100.0%', result.stderr)
        self.assertNotIn('Write retained markers:', result.stderr)

    def test_random_string_oracle_all_bit_boundaries_and_threads(self):
        rng = random.Random(7401)
        sequences = [''.join(rng.choices("ACGT", k=180)) for _ in range(5)]
        sequences[2] = sequences[2][:90] + "NNrySW" + sequences[2][90:]
        sequences[3] = sequences[3].lower()
        target = ''.join(row(f"t{i}", seq, "9000001") for i, seq in enumerate(sequences))
        background_seqs = [reverse_complement(sequences[0][20:130]), sequences[1][30:120], "N"*80]
        background = target + ''.join(row(f"b{i}", seq, "9606") for i, seq in enumerate(background_seqs))
        for k in (1, 2, 7, 31, 32, 33, 42, 63, 64):
            with self.subTest(k=k):
                expected = set().union(*(kmers(s, k) for s in sequences)) - set().union(*(kmers(s, k) for s in background_seqs))
                first, _, _ = self.execute(target, background, k=k)
                parallel, _, _ = self.execute(target, background, k=k, threads=4)
                self.assertEqual(first, parallel)
                self.assertEqual({r["sequence"] for r in first}, expected)
                for r in first:
                    source = sequences[int(r["source_accession"][1:])].upper()
                    start = int(r["start_1based"])-1
                    window = source[start:start+k]
                    self.assertEqual(r["sequence"], reverse_complement(window) if r["strand"] == '-' else window)

    def test_42mer_high_bits_are_not_truncated(self):
        suffix = "CGTAGCTACGATCGTAGCATCGTAGCTACGTT"
        a, b = "AAAAAAAAAA"+suffix, "AAAAAAAACA"+suffix
        records, _, _ = self.execute(row("a", a)+row("b", b), row("bg", a, "9606"))
        self.assertEqual([r["sequence"] for r in records], [canonical(b)])

    def test_mixed_unknown_and_merged_deflines_are_background(self):
        rng = random.Random(7)
        seqs = [''.join(rng.choices("ACGT", k=42)) for _ in range(10)]
        target = ''.join(row(f"t{i}", s) for i, s in enumerate(seqs))
        background = target
        background += row("mixed", seqs[0], "10226", "10227;9606")
        background += row("unknown", seqs[1], "N/A")
        background += row("unknown0", seqs[2], "0")
        background += row("other_defline_same_sequence", reverse_complement(seqs[3]), "9606")
        background += row("ancestor", seqs[4], "33208", "10226")
        background += row("unmapped", seqs[5], "4000000000")
        background += row("mixed_principal", seqs[6], "10226;9606")
        background += row("mixed_unknown_leaf", seqs[7], "10226", "10227;0")
        background += row("all_target", seqs[8], "10226", "10227;9000001")
        background += row("all_target_0_leaves", seqs[9], "10227", "0")
        records, _, _ = self.execute(target, background, threads=30)
        self.assertEqual({r["sequence"] for r in records}, {canonical(s) for s in seqs[8:]})

    def test_no_crossing_record_or_ambiguous_base_boundaries(self):
        sequence = "ACGTCGATGCATCGTAGCATGCTACGTACGTCGATCGTAGCA"
        background = row("left", sequence[:21], "9606") + row("right", sequence[21:], "9606")
        background += row("ambiguous", sequence[:21]+"N"+sequence[21:], "9606")
        records, _, _ = self.execute(row("target", sequence), background)
        self.assertEqual([r["sequence"] for r in records], [canonical(sequence)])

    def test_chunk_overlap_at_one_megabyte(self):
        sequence = "ACGTCGATGCATCGTAGCATGCTACGTACGTCGATCGTAGCA"
        background = "N"*(1024*1024-20) + reverse_complement(sequence) + "N"*100
        records, _, _ = self.execute(row("target", sequence), row("bg", background, "9606"), threads=4)
        self.assertEqual(records, [])

    def test_all_shared_is_successful_empty_result(self):
        sequence = "ACGT"*20
        records, manifest, output = self.execute(row("target", sequence), row("bg", sequence, "9606"))
        self.assertEqual(records, [])
        self.assertEqual((output/"markers.fasta").read_text(), "")
        self.assertEqual(manifest["removed_kmers"], manifest["candidate_kmers"])

    def test_failed_exports_are_not_accepted(self):
        target, background = row("target", "ACGT"*20), row("bg", "A"*80, "9606")
        for config in ({"target_exit": 8}, {"background_exit": 9}, {"changed": True}):
            with self.subTest(config=config):
                result, output = self.execute(target, background, config=config, success=False)
                self.assertIn("ERROR:", result.stderr)
                self.assertIn("status\tfailed", (output/"manifest.tsv").read_text())
                self.assertFalse((output/"markers.fasta").exists())

    def test_empty_truncated_malformed_or_uninformative_exports(self):
        target = row("target", "ACGT"*20)
        for background in ("", "broken\n", row("b", "A"*80, "9606").rstrip("\n"), target):
            with self.subTest(background=background[:20]):
                self.execute(target, background, success=False)
        self.execute(row("short", "ACT"), row("b", "A"*80, "9606"), success=False)

    def test_taxonomy_validation(self):
        original = (self.root/"nodes.dmp").read_text()
        malformed = [original.replace("10227\t|\t10226", "10227\t|\t9000001"),
                     original.replace("10226\t|\t33208", "10226\t|\t444444"),
                     original + "10226\t|\t33208\t|\n"]
        for nodes in malformed:
            (self.root/"nodes.dmp").write_text(nodes)
            self.execute(row("t", "ACGT"*20), row("b", "A"*80, "9606"), success=False)
        (self.root/"nodes.dmp").write_text(original)
        (self.root/"names.dmp").write_text("1\t|\troot\t|\t\t|\tscientific name\t|\n")
        self.execute(row("t", "ACGT"*20), row("b", "A"*80, "9606"), success=False)

    def test_cli_errors_and_existing_output(self):
        for args in (["-k", "65"], ["-k", "0"], ["-t", "0"], ["-k", "x"],
                     ["-k", "42", "-k", "42"], ["-taxid", "-1"],
                     ['-progress-seconds','0'], ['-progress-seconds','3601'],
                     ['-progress-seconds','-1'], ['-progress-seconds','1.5']):
            result = subprocess.run([str(BINARY)]+args, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
        output = self.root/"output 1"
        output.mkdir()
        sentinel = output/"user_file"
        sentinel.write_text("keep")
        self.execute(row("t", "ACGT"*20), row("b", "A"*80, "9606"), success=False)
        self.assertEqual(sentinel.read_text(), "keep")

    @unittest.skipUnless(SQLITE3_BINARY, 'sqlite3 CLI required for BLAST taxonomy tests')
    def test_blast_taxonomy_equivalence_and_provenance(self):
        taxonomy = write_blast_taxonomy(self.root/"taxonomy ' ; $(touch TAX_INJECTED)", self.root)
        rng = random.Random(197)
        sequences = [''.join(rng.choices('ACGT', k=80)) for _ in range(6)]
        target = ''.join(row('t'+str(i), seq, '10226' if i % 2 else '9000001')
                         for i, seq in enumerate(sequences))
        # Unknown/mixed/ancestral principals remain background in either mode.
        background = target + row('mixed', sequences[0], '10226', '10227;9606')
        background += row('unknown', sequences[1], '0')
        background += row('ancestor', sequences[2], '33208', '10227')
        background += row('reverse', reverse_complement(sequences[3]), '9606')
        expected = set().union(*(kmers(s,42) for s in sequences))
        expected -= set().union(*(kmers(s,42) for s in sequences[:4]))
        dumps, _, dump_output = self.execute(target, background)
        for workers in (1, 30):
            records, manifest, output = self.execute(target, background, blasttax=taxonomy, threads=workers)
            self.assertEqual(records, dumps)
            self.assertEqual({r['sequence'] for r in records}, expected)
            self.assertEqual(manifest['target_name'], 'Placozoa')
            self.assertEqual(manifest['taxonomy_nodes'], '7')
            self.assertIn('sqlite3_version', manifest)
            for key, name in [('taxonomy_sqlite','taxonomy4blast.sqlite3'),
                              ('taxdb_bti','taxdb.bti'), ('taxdb_btd','taxdb.btd')]:
                self.assertEqual(manifest[key+'_file'], str(taxonomy/name))
                self.assertEqual(int(manifest[key+'_bytes']), (taxonomy/name).stat().st_size)
                self.assertRegex(manifest[key+'_fnv1a64'], r'^[0-9a-f]{16}$')
            self.assertEqual((output/'target_taxids.txt').read_bytes(), (dump_output/'target_taxids.txt').read_bytes())
            self.assertTrue((output/'taxonomy.log').exists())
            self.assertIn('-readonly', (output/'commands.log').read_text())

    @unittest.skipUnless(SQLITE3_BINARY, 'sqlite3 CLI required for BLAST taxonomy tests')
    def test_blast_taxonomy_root_leaf_and_name_index_boundaries(self):
        taxonomy = write_blast_taxonomy(self.root/'taxonomy', self.root)
        for taxid, ids in [(1, [1,2,9606,10226,10227,33208,9000001]), (9000001,[9000001])]:
            with self.subTest(taxid=taxid):
                records, manifest, _ = self.execute(row('target','ACGT'*20,'9000001'),
                    row('unknown','A'*80,'0'), blasttax=taxonomy, target_taxid=taxid,
                    config={'expected_taxids':list(map(str,ids))})
                self.assertTrue(records)
                self.assertEqual(manifest['target_name'], 'taxon '+str(taxid))

    def test_taxonomy_mode_cli_validation(self):
        base = [str(BINARY), '-taxid','10226','-blastdb','unused']
        for args in ([], ['-taxnodes','nodes.dmp'], ['-taxnames','names.dmp'],
                     ['-blasttax','taxonomy','-taxnodes','nodes.dmp'],
                     ['-blasttax','taxonomy','-taxnames','names.dmp'],
                     ['-taxnodes','nodes.dmp','-taxnames','names.dmp','-sqlite3','sqlite3']):
            with self.subTest(args=args):
                result = subprocess.run(base+args, cwd=self.root, capture_output=True, text=True)
                self.assertNotEqual(result.returncode,0)
                self.assertFalse((self.root/'clademers_10226_k42').exists())

    @unittest.skipUnless(SQLITE3_BINARY, 'sqlite3 CLI required for BLAST taxonomy tests')
    def test_blast_taxonomy_missing_corrupt_and_malformed_files(self):
        cases = ['taxonomy4blast.sqlite3','taxdb.bti','taxdb.btd','magic','count','offset','order',
                 'name','empty_name','sqlite','wal','journal']
        for case in cases:
            with self.subTest(case=case):
                directory = write_blast_taxonomy(self.root/case, self.root)
                if case in ('taxonomy4blast.sqlite3','taxdb.bti','taxdb.btd'):
                    (directory/case).unlink()
                elif case in ('wal','journal'):
                    (directory/('taxonomy4blast.sqlite3-'+case)).write_text('active transaction')
                elif case == 'sqlite':
                    (directory/'taxonomy4blast.sqlite3').write_text('not a database')
                elif case in ('name','empty_name'):
                    content = (directory/'taxdb.btd').read_bytes()
                    content = content.replace(b'Placozoa', b'\0lacozoa' if case=='name' else b'        ')
                    (directory/'taxdb.btd').write_bytes(content)
                else:
                    content = bytearray((directory/'taxdb.bti').read_bytes())
                    offsets = {'magic':0, 'count':4, 'offset':28, 'order':32}
                    values = {'magic':1, 'count':999, 'offset':2**32-1, 'order':1}
                    struct.pack_into('>I',content,offsets[case],values[case])
                    (directory/'taxdb.bti').write_bytes(content)
                result, output = self.execute(row('t','ACGT'*20),row('b','A'*80,'9606'),
                                               blasttax=directory,success=False)
                self.assertIn('ERROR:', result.stderr)
                self.assertIn('status\tfailed', (output/'manifest.tsv').read_text())

    @unittest.skipUnless(SQLITE3_BINARY, 'sqlite3 CLI required for BLAST taxonomy tests')
    def test_blast_taxonomy_missing_internal_name_uses_taxid_label(self):
        directory = write_blast_taxonomy(self.root/'taxonomy', self.root)
        index = bytearray((directory/'taxdb.bti').read_bytes())
        struct.pack_into('>I',index,48,10225)  # Keep the index sorted, omit 10226.
        (directory/'taxdb.bti').write_bytes(index)
        records, manifest, _ = self.execute(row('t','ACGT'*20),row('b','A'*80,'9606'),blasttax=directory)
        self.assertTrue(records)
        self.assertEqual(manifest['target_name'],'taxid 10226')
        self.assertEqual(manifest['target_name_source'],'taxid_only')
        # A missing tree node is still an error, even when a numeric label exists.
        self.execute(row('t','ACGT'*20),row('b','A'*80,'9606'),blasttax=directory,
                     target_taxid=444444,success=False)

    @unittest.skipUnless(SQLITE3_BINARY, 'sqlite3 CLI required for BLAST taxonomy tests')
    def test_blast_taxonomy_invalid_trees_and_schema(self):
        edits = ["DROP TABLE TaxidInfo", "ALTER TABLE TaxidInfo RENAME COLUMN parent TO wrong",
                 "DELETE FROM TaxidInfo", "DELETE FROM TaxidInfo WHERE taxid=10226",
                 "UPDATE TaxidInfo SET parent=9000001 WHERE taxid=10227",
                 "UPDATE TaxidInfo SET parent=444444 WHERE taxid=10226",
                 "UPDATE TaxidInfo SET parent=NULL WHERE taxid=10226",
                 "UPDATE TaxidInfo SET parent='invalid' WHERE taxid=10226",
                 "UPDATE TaxidInfo SET parent=2 WHERE taxid=2",
                 "CREATE TABLE duplicate AS SELECT * FROM TaxidInfo; DROP TABLE TaxidInfo; "
                 "ALTER TABLE duplicate RENAME TO TaxidInfo; INSERT INTO TaxidInfo VALUES(10226,33208)"]
        for i, sql in enumerate(edits):
            with self.subTest(sql=sql):
                directory = write_blast_taxonomy(self.root/('bad-tree-'+str(i)), self.root)
                with sqlite3.connect(directory/'taxonomy4blast.sqlite3') as connection:
                    connection.executescript(sql)
                self.execute(row('t','ACGT'*20),row('b','A'*80,'9606'),blasttax=directory,success=False)

    def test_blast_taxonomy_failed_and_truncated_sqlite_export(self):
        directory = write_blast_taxonomy(self.root/'taxonomy', self.root)
        nodes = ''.join(str(t)+'\t'+str(p)+'\n' for t,p in
                        [(1,1),(2,1),(33208,1),(10226,33208),(10227,10226),(9000001,10227),(9606,33208)])
        wrapper = self.root/'sqlite reader.py'
        for exported, code in [('ok\n'+nodes,9), ('ok\n'+nodes.rstrip('\n'),0), ('ok\nbroken\n',0)]:
            wrapper.write_text('#!'+sys.executable+'\nimport sys\n'
                'if "-version" in sys.argv: print("synthetic sqlite reader")\n'
                'else:\n    sys.stdout.write('+repr(exported)+')\n    sys.exit('+str(code)+')\n')
            wrapper.chmod(0o755)
            self.execute(row('t','ACGT'*20),row('b','A'*80,'9606'),blasttax=directory,
                         sqlite3_reader=str(wrapper),success=False)
        result, _ = self.execute(row('t','ACGT'*20),row('b','A'*80,'9606'),blasttax=directory,
                                sqlite3_reader=str(self.root/'missing-sqlite3'),success=False)
        self.assertIn('Cannot start', result.stderr)

    @unittest.skipUnless(SQLITE3_BINARY, 'sqlite3 CLI required for BLAST taxonomy tests')
    def test_taxonomy_changes_during_run_are_rejected(self):
        directory = write_blast_taxonomy(self.root/'taxonomy', self.root)
        for blasttax, path in [(None,self.root/'names.dmp'), (directory,directory/'taxdb.btd')]:
            result, output = self.execute(row('t','ACGT'*20),row('b','A'*80,'9606'),
                blasttax=blasttax, config={'change_taxonomy_file':str(path)}, success=False)
            self.assertIn('Taxonomy file changed',result.stderr)
            self.assertFalse((output/'markers.fasta').exists())

    @unittest.skipUnless(os.environ.get('BLAST_TAXDB') and os.environ.get('BLAST_BIN') and SQLITE3_BINARY,
                         'set BLAST_TAXDB and BLAST_BIN, with sqlite3 CLI, for an NCBI taxonomy snapshot')
    def test_real_ncbi_blast_taxonomy(self):
        directory = Path(os.environ['BLAST_TAXDB']).resolve()
        bindir = Path(os.environ['BLAST_BIN'])
        # A separate recursive SQL oracle defines the expected descendants.
        with sqlite3.connect((directory/'taxonomy4blast.sqlite3').as_uri()+'?mode=ro', uri=True) as con:
            ids = [r[0] for r in con.execute('WITH RECURSIVE clade(taxid) AS ('
                'SELECT taxid FROM TaxidInfo WHERE taxid=10226 UNION '
                'SELECT t.taxid FROM TaxidInfo t JOIN clade c ON t.parent=c.taxid) '
                'SELECT taxid FROM clade ORDER BY taxid')]
            # Derive a small .dmp tree containing this entire clade, the demo's
            # background taxa, and all their ancestors. No taxdump download.
            parents = {}
            for start in ids + [9606,2]:
                taxid = start
                while taxid not in parents:
                    parent = con.execute('SELECT parent FROM TaxidInfo WHERE taxid=?',(taxid,)).fetchone()[0]
                    parents[taxid] = parent
                    taxid = parent
        self.assertIn(10227, ids)
        (self.root/'nodes.dmp').write_text(''.join(f'{t}\t|\t{p}\t|\tno rank\t|\n'
                                                  for t,p in sorted(parents.items())))
        (self.root/'names.dmp').write_text('10226\t|\tPlacozoa\t|\t\t|\tscientific name\t|\n')
        database = self.root/'demo_db'
        made = subprocess.run([str(bindir/'makeblastdb'),'-in',str(ROOT/'examples/references.fasta'),
            '-dbtype','nucl','-parse_seqids','-taxid_map',str(ROOT/'examples/taxid_map.tsv'),
            '-blastdb_version','5','-out',str(database)],capture_output=True,text=True)
        self.assertEqual(made.returncode,0,made.stderr)
        config = {'expected_taxids':list(map(str,ids))}
        dumps, _, _ = self.execute('','',database=database,blast=bindir/'blastdbcmd',config=config)
        records, manifest, _ = self.execute('','',database=database,blast=bindir/'blastdbcmd',
                                             blasttax=directory,config=config,threads=30)
        self.assertEqual(records,dumps)
        self.assertIn(manifest['target_name'],('Placozoa','taxid 10226'))
        self.assertEqual(manifest['target_name_source'],
                         'taxid_only' if manifest['target_name']=='taxid 10226' else 'taxdb')
        self.assertEqual(manifest['candidate_kmers'],'2')
        self.assertEqual(manifest['retained_kmers'],'1')
        for source in (None, directory):
            fasta_records, _, _ = self.execute_fasta((ROOT/'examples/references.fasta').read_text(),
                (ROOT/'examples/taxid_map.tsv').read_text(), blasttax=source)
            self.assertEqual(common_provenance(fasta_records), common_provenance(records))
        # Exercise a real name record too: some releases omit Placozoa's name.
        with sqlite3.connect((directory/'taxonomy4blast.sqlite3').as_uri()+'?mode=ro', uri=True) as con:
            human_ids = [r[0] for r in con.execute('WITH RECURSIVE clade(taxid) AS ('
                'SELECT taxid FROM TaxidInfo WHERE taxid=9606 UNION '
                'SELECT t.taxid FROM TaxidInfo t JOIN clade c ON t.parent=c.taxid) '
                'SELECT taxid FROM clade ORDER BY taxid')]
        human_records, human_manifest, _ = self.execute('','',database=database,blast=bindir/'blastdbcmd',
            blasttax=directory,target_taxid=9606,config={'expected_taxids':list(map(str,human_ids))})
        self.assertEqual(human_manifest['target_name'],'Homo sapiens')
        self.assertEqual(human_manifest['target_name_source'],'taxdb')
        # The sole human sequence in this fixture is also in Placozoa.
        self.assertEqual(human_records,[])
        self.assertEqual(human_manifest['candidate_kmers'],'1')
        self.assertEqual(human_manifest['retained_kmers'],'0')

    @unittest.skipUnless(os.environ.get("BLAST_BIN"), "set BLAST_BIN for genuine BLAST+ integration")
    def test_real_blast_version5_database(self):
        bindir = Path(os.environ["BLAST_BIN"])
        rng = random.Random(987)
        sequences = [''.join(rng.choices("ACGT", k=150)) for _ in range(5)]
        # Two taxonomically distinct records have the exact same sequence;
        # regardless of BLAST's internal storage, this sequence must be excluded.
        entries = [("target_a", sequences[0], 10226), ("target_b", sequences[1], 10227),
                   ("target_c", sequences[2], 9000001), ("bg_a", sequences[0], 9606),
                   ("bg_b", reverse_complement(sequences[1][20:110]), 9606),
                   ("bg_c", sequences[3], 2), ("unknown", sequences[4], 0)]
        fasta, taxmap = self.root/"refs.fasta", self.root/"taxmap.tsv"
        fasta.write_text(''.join(f">{name}\n{seq}\n" for name, seq, _ in entries))
        taxmap.write_text(''.join(f"{name}\t{tax}\n" for name, _, tax in entries))
        # makeblastdb 2.17 refuses spaces in -out even with inner quotes.
        # Build in a plain directory, then move it to test blastdbcmd's documented
        # quoted database-path support independently of makeblastdb's restriction.
        staging = self.root/"build"
        staging.mkdir()
        database = staging/"db"
        made = subprocess.run([str(bindir/"makeblastdb"), "-in", str(fasta), "-dbtype", "nucl",
                               "-parse_seqids", "-taxid_map", str(taxmap), "-blastdb_version", "5",
                               "-out", str(database)], capture_output=True, text=True)
        self.assertEqual(made.returncode, 0, made.stderr)
        staging.rename(self.root/"real database")
        database = self.root/"real database"/"db"
        expected = set().union(*(kmers(s, 42) for _, s, t in entries if t in (10226, 10227, 9000001)))
        expected -= set().union(*(kmers(s, 42) for _, s, t in entries if t not in (10226, 10227, 9000001)))
        serial, _, _ = self.execute("", "", database=database, blast=bindir/"blastdbcmd")
        self.assertIn('BLAST background scan: 0.0%', self.last_result.stderr)
        self.assertIn('BLAST background scan: 100.0%', self.last_result.stderr)
        parallel, _, _ = self.execute("", "", database=database, blast=bindir/"blastdbcmd", threads=30)
        self.assertEqual({r["sequence"] for r in serial}, expected)
        self.assertEqual(serial, parallel)

        if GZIP_BINARY:
            for nodes,names in ((True,False),(False,True),(True,True)):
                records, _, _ = self.execute('', '',database=database,blast=bindir/'blastdbcmd',
                    nodes_compressed=nodes,names_compressed=names,threads=4)
                self.assertEqual(records,serial)

        headers = ''.join(f'>{name} taxid={tax}\n{seq}\n' for name,seq,tax in entries)
        for mapped in (False,True):
            records, _, _ = self.execute_fasta(fasta.read_text() if mapped else headers,
                                               taxmap.read_text() if mapped else None)
            self.assertEqual(common_provenance(records), common_provenance(serial))

        if SQLITE3_BINARY:
            taxonomy = write_blast_taxonomy(self.root/'blast taxonomy', self.root)
            for workers in (1,30):
                records, _, _ = self.execute('', '', database=database, blast=bindir/'blastdbcmd',
                                              blasttax=taxonomy, threads=workers)
                self.assertEqual(records,serial)
                for mapped in (False,True):
                    records, _, _ = self.execute_fasta(fasta.read_text() if mapped else headers,
                        taxmap.read_text() if mapped else None, blasttax=taxonomy, threads=workers)
                    self.assertEqual(common_provenance(records), common_provenance(serial))


if __name__ == "__main__":
    if not BINARY.is_file():
        sys.exit("Build CladeMers first: make")
    unittest.main(verbosity=2)
