#!/usr/bin/env python3
"""Direct BLAST/plain/indexed-gzip restart tests; no reference downloads."""
import gzip
import json
import os
from pathlib import Path
import random
import re
import shutil
import signal
import struct
import subprocess
import sys
import time
import unittest
import zlib

import test_checkpoint as checkpoints
from test_clademers import BINARY, ROOT, fnv1a64


class DirectResumeTests(unittest.TestCase):
    setUp = checkpoints.CheckpointTests.setUp
    tearDown = checkpoints.CheckpointTests.tearDown
    configure = checkpoints.CheckpointTests.configure
    run_ok = checkpoints.CheckpointTests.run_ok
    run_bad = checkpoints.CheckpointTests.run_bad
    launch = checkpoints.CheckpointTests.launch
    wait_log = checkpoints.CheckpointTests.wait_log
    interrupt = checkpoints.CheckpointTests.interrupt
    manifest = checkpoints.CheckpointTests.manifest
    compare = checkpoints.CheckpointTests.compare
    checkpoint_digest = checkpoints.CheckpointTests.checkpoint_digest

    def fasta_args(self, path, output):
        return ['-fasta',str(path),'-taxnodes',str(self.root/'nodes.dmp'),
                '-taxnames',str(self.root/'names.dmp'),'-taxid','10226','-out',str(output),
                '-checkpoint-seconds','1','-t','4']

    def test_indexed_gzip_headers_members_boundaries_and_corruption(self):
        content=(ROOT/'examples/references_taxids.fasta').read_bytes().replace(b'\n',b'\r\n').rstrip(b'\r\n')
        plain=self.root/'plain.fa'; plain.write_bytes(content)
        baseline=self.root/'baseline'; self.run_ok(self.fasta_args(plain,baseline))
        header=b'\x1f\x8b\x08\x1e'+b'\0'*4+b'\0\xff'+struct.pack('<H',3)+b'abc'+b'references.fa\0comment\0'
        header+=struct.pack('<H',zlib.crc32(header)&65535)
        def member(data,level=6,strategy=zlib.Z_DEFAULT_STRATEGY):
            compressor=zlib.compressobj(level,zlib.DEFLATED,-15,8,strategy)
            return header+compressor.compress(data)+compressor.flush()+struct.pack('<II',zlib.crc32(data),len(data))
        valid=[gzip.compress(content,mtime=0),member(content),member(content,0),member(content,6,zlib.Z_FIXED),
               member(content[:31])+member(b'')+member(content[31:]),gzip.compress(content)+b'\0'*19]
        for i,data in enumerate(valid):
            with self.subTest(valid=i):
                path=self.root/(str(i)+'.payload'); path.write_bytes(data)
                output=self.root/('valid-'+str(i)); self.run_ok(self.fasta_args(path,output))
                self.assertEqual((baseline/'markers.tsv').read_bytes(),(output/'markers.tsv').read_bytes())
                manifest=self.manifest(output)
                self.assertEqual(manifest['fasta_reader'],'indexed_zlib')
                self.assertEqual(manifest['fasta_fnv1a64'],fnv1a64(content))
                self.assertFalse((output/'gzip.log').exists())
        good=bytearray(member(content))
        bad_crc=bytearray(good); bad_crc[-8]^=1
        bad_size=bytearray(good); bad_size[-4]^=1
        bad_header=bytearray(good); bad_header[len(header)-1]^=1
        bad_flags=bytearray(good); bad_flags[3]|=0x80
        invalid=[bytes(good[:2]),bytes(good[:12]),bytes(good[:-1]),bytes(good[:-8]),bytes(bad_crc),bytes(bad_size),
                 bytes(bad_header),bytes(bad_flags),bytes(good)+b'garbage',bytes(good)+b'\0x',bytes(good)+member(b'')[:-2]]
        for i,data in enumerate(invalid):
            with self.subTest(invalid=i):
                path=self.root/('bad-'+str(i)+'.gz'); path.write_bytes(data)
                output=self.root/('bad-result-'+str(i))
                self.run_bad(self.fasta_args(path,output),'gzip')
                self.assertFalse((output/'SUCCESS').exists())

    def large_fasta(self, path, phase, compressed, concatenated=False):
        # Low candidate cardinality, enough work for a timed periodic checkpoint.
        # Stream generation so the fixture doesn't require a second large RAM copy.
        rng=random.Random(782)
        motif=''.join(rng.choices('ACGT',k=256)).encode()
        sequence=motif*4096 # 1 MiB record, many gzip DEFLATE block boundaries overall.
        with path.open('wb') as raw:
            def member(): return gzip.GzipFile(fileobj=raw,mode='wb',compresslevel=1)
            out=member() if compressed else raw
            try:
                if phase=='background scan': out.write(b'>target taxid=10226\n'+motif*2+b'N'+b'C'*42+b'\n')
                for i in range(512 if phase=='background scan' else 128):
                    taxid=b'10226' if phase=='target loading' else b'9606'
                    out.write(b'>s'+str(i).encode()+b' taxid='+taxid+b'\n'+sequence[:123])
                    if compressed and concatenated and i%4==0:
                        out.close(); out=member(); out.close(); out=member() # Empty member, inside a sequence line.
                    out.write(sequence[123:]+b'\n')
                if phase=='target loading': out.write(b'>background taxid=9606\n'+b'C'*42)
            finally:
                if compressed: out.close()

    def direct_fasta_restart(self, phase, compressed, other=None, concatenated=False):
        path=self.root/('large-'+phase.replace(' ','-')+('.fa.gz' if compressed else '.fa'))
        self.large_fasta(path,phase,compressed,concatenated)
        baseline=self.root/('baseline-'+phase.replace(' ','-'))
        self.run_ok(self.fasta_args(path,baseline))
        output=self.root/('paused-'+phase.replace(' ','-'))
        process,log=self.launch(self.fasta_args(path,output))
        self.wait_log(process,log,'Saved '+phase+r' after [1-9][0-9]* records',regex=True)
        self.interrupt(process,log,signal.SIGKILL)
        result=self.run_ok(['-resume',str(output),'-t','1'],binary=other or BINARY)
        self.assertNotIn('Replay saved input prefix',result.stderr)
        if compressed:
            match=re.search(r'Indexed gzip seek: compressed byte (\d+); replayed (\d+) content bytes',result.stderr)
            self.assertIsNotNone(match,result.stderr)
            self.assertGreater(int(match.group(1)),10)
            # Standard compressor blocks need only a small local replay, not the saved prefix.
            self.assertLess(int(match.group(2)),8*1024*1024)
        else:
            self.assertRegex(result.stderr,r'Plain FASTA seek: byte [1-9][0-9]+')
        self.compare(baseline,output)
        self.assertEqual(self.manifest(baseline)['fasta_fnv1a64'],self.manifest(output)['fasta_fnv1a64'])
        self.assertEqual(self.manifest(baseline)['fasta_records'],self.manifest(output)['fasta_records'])
        return path,output

    def test_plain_fasta_seeks_to_saved_record(self):
        self.direct_fasta_restart('target loading',False)

    def test_indexed_gzip_target_and_background_seek(self):
        for phase in ('target loading','background scan'):
            with self.subTest(phase=phase): self.direct_fasta_restart(phase,True)

    def test_indexed_gzip_seek_after_concatenated_and_empty_members(self):
        self.direct_fasta_restart('target loading',True,concatenated=True)

    def test_changed_file_even_with_restored_size_and_mtime_is_rejected(self):
        path=self.root/'input.fa'; path.write_bytes((ROOT/'examples/references_taxids.fasta').read_bytes())
        baseline=self.root/'baseline'; self.run_ok(self.fasta_args(path,baseline))
        output=self.root/'paused'; shutil.copytree(baseline,output); (output/'SUCCESS').unlink()
        before=path.stat(); data=path.read_bytes(); path.write_bytes(data.replace(b'A',b'C',1))
        os.utime(path,ns=(before.st_atime_ns,before.st_mtime_ns))
        self.run_bad(['-resume',str(output)],'Checkpoint inputs changed')
        self.assertFalse((output/'SUCCESS').exists())

    def real_db(self):
        bindir=Path(os.environ['BLAST_BIN']); real=bindir/'blastdbcmd'
        directory=self.root/'real-build'; directory.mkdir()
        rng=random.Random(873)
        for n in range(2):
            sub=directory/('vol'+str(n)); sub.mkdir()
            entries=[('seq'+str(n)+'_'+str(i),''.join(rng.choices('ACGT',k=100)),10226 if i%2==0 else 9606) for i in range(8)]
            fasta=sub/'references.fa'; mapping=sub/'taxids.tsv'
            fasta.write_text(''.join('>'+name+'\n'+seq+'\n' for name,seq,_ in entries))
            mapping.write_text(''.join(name+'\t'+str(tax)+'\n' for name,_,tax in entries))
            made=subprocess.run([str(bindir/'makeblastdb'),'-in',str(fasta),'-dbtype','nucl','-parse_seqids',
                '-taxid_map',str(mapping),'-out',str(sub/'db')],capture_output=True,text=True)
            self.assertEqual(made.returncode,0,made.stderr)
        directory.rename(self.root/"real database ' ; $(touch INJECTED)")
        directory=self.root/"real database ' ; $(touch INJECTED)"
        (directory/'all.nal').write_text('TITLE two volumes\nDBLIST "vol0/db" "vol1/db"\n')
        (directory/'selected.nal').write_text('TITLE original filter retained\nDBLIST "all"\nFIRST_OID 2\nLAST_OID 15\n')
        return real,directory/'selected'

    @unittest.skipUnless(os.environ.get('BLAST_BIN'),'real BLAST+ required')
    def test_blast_seeks_global_oids_across_nested_aliases_and_deflines(self):
        real,database=self.real_db()
        audit=self.root/'exported-oids.txt'; delay=self.root/'delay'; delay.write_text('0')
        wrapper=self.root/'observed-blast.py'
        wrapper.write_text('#!'+sys.executable+'\n'+
            'import pathlib,subprocess,sys,time\n'+
            'args=sys.argv[1:]; slow=float(pathlib.Path('+repr(str(delay))+').read_text())\n'+
            'p=subprocess.Popen(['+repr(str(real))+']+args,stdout=subprocess.PIPE)\n'+
            'for line in p.stdout:\n'+
            ' copies=2 if "-outfmt" in args else 1\n'+
            ' for copy in range(copies):\n'+
            '  if "-outfmt" in args:\n'+
            '   with open('+repr(str(audit))+',"a") as f: f.write(line.split(b"\\t",1)[0].decode()+"\\n")\n'+
            '  sys.stdout.buffer.write(line); sys.stdout.buffer.flush()\n'+
            '  if "-outfmt" in args: time.sleep(slow)\n'+
            'sys.exit(p.wait())\n')
        wrapper.chmod(0o755)
        common=['-blastdb',str(database),'-blastdbcmd',str(wrapper),'-taxnodes',str(self.root/'nodes.dmp'),
                '-taxnames',str(self.root/'names.dmp'),'-taxid','10226','-checkpoint-seconds','1','-t','4']
        baseline=self.root/'baseline'; self.run_ok(common+['-out',str(baseline)])
        for phase in ('target loading','background scan'):
            with self.subTest(phase=phase):
                delay.write_text('0.28'); audit.write_text('')
                output=self.root/phase.replace(' ','-')
                process,log=self.launch(common+['-out',str(output)])
                self.wait_log(process,log,'Saved '+phase+r' after [1-9][0-9]* records',regex=True)
                self.interrupt(process,log)
                delay.write_text('0'); audit.write_text('')
                result=self.run_ok(['-resume',str(output),'-t','1'])
                match=re.search(r'Direct BLAST seek to OID (\d+); validating (\d+) saved',result.stderr)
                self.assertIsNotNone(match,result.stderr)
                start=int(match.group(1)); self.assertGreater(start,1)
                exported=list(map(int,audit.read_text().splitlines()))
                self.assertEqual(exported[0],start)
                # Target resume is followed by an intentional full background pass.
                if phase=='background scan': self.assertTrue(all(oid>=start for oid in exported))
                self.assertIn('Saved OID group verified',result.stderr)
                self.assertNotIn('Replay saved input prefix',result.stderr)
                self.compare(baseline,output)
                self.assertFalse((self.root/'INJECTED').exists())

    @unittest.skipUnless(os.environ.get('BLAST_BIN'),'real BLAST+ required')
    def test_changed_blast_volume_and_alias_rejected_before_export(self):
        real,database=self.real_db()
        common=['-blastdb',str(database),'-blastdbcmd',str(real),'-taxnodes',str(self.root/'nodes.dmp'),
                '-taxnames',str(self.root/'names.dmp'),'-taxid','10226']
        for i,path in enumerate((database.with_suffix('.nal'),database.parent/'vol0/db.nsq')):
            baseline=self.root/('baseline-'+str(i)); self.run_ok(common+['-out',str(baseline)])
            output=self.root/('changed-'+str(i)); shutil.copytree(baseline,output); (output/'SUCCESS').unlink()
            # A harmless timestamp change alone must invalidate a direct-seek snapshot.
            os.utime(path,None)
            self.run_bad(['-resume',str(output)],'Checkpoint inputs changed')

    @unittest.skipUnless(os.environ.get('CLADEMERS_RESUME_BINARY'),'second-architecture binary required')
    def test_gzip_index_checkpoint_on_second_architecture(self):
        self.direct_fasta_restart('background scan',True,Path(os.environ['CLADEMERS_RESUME_BINARY']))


if __name__=='__main__': unittest.main(verbosity=2)
