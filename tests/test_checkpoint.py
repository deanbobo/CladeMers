#!/usr/bin/env python3
"""Restart tests use real process signals and compare complete output/provenance."""
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import random
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import unittest

from test_clademers import BINARY, FAKE_BLAST, GZIP_BINARY, SQLITE3_BINARY, ROOT, row, write_taxonomy, write_blast_taxonomy


class CheckpointTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='clademers-restart-')
        self.root = Path(self.temp.name)
        write_taxonomy(self.root)
        self.fake = self.root/'blast-reader.py'
        self.fake.write_text('#!'+sys.executable+'\n'+FAKE_BLAST)
        self.fake.chmod(0o755)
        self.db = self.root/'db'
        self.db.mkdir()
        rng = random.Random(439)
        self.sequences = [''.join(rng.choices('ACGT',k=100)) for _ in range(8)]
        self.target = ''.join(row('t'+str(i),s) for i,s in enumerate(self.sequences))
        self.background = self.target+''.join(row('b'+str(i),s,'9606') for i,s in enumerate(self.sequences[:4]))
        (self.db/'target.tsv').write_text(self.target)
        (self.db/'all.tsv').write_text(self.background)
        self.configure()
        self.children = []
        self.serial = 0

    def tearDown(self):
        for process in self.children:
            if process.poll() is None:
                os.killpg(process.pid,signal.SIGKILL)
                process.wait(timeout=10)
        self.temp.cleanup()

    def configure(self, **values):
        (self.db/'config.json').write_text(json.dumps(values))

    def args(self, output, k=42):
        return ['-resume-mode','replay','-blastdb',str(self.db),'-blastdbcmd',str(self.fake),
                '-taxnodes',str(self.root/'nodes.dmp'),'-taxnames',str(self.root/'names.dmp'),
                '-taxid','10226','-k',str(k),'-out',str(output),'-checkpoint-seconds','1']

    def run_ok(self, args, binary=BINARY):
        result = subprocess.run([str(binary)]+args,cwd=self.root,capture_output=True,text=True,timeout=60)
        self.assertEqual(result.returncode,0,result.stderr)
        return result

    def run_bad(self, args, message=None):
        result = subprocess.run([str(BINARY)]+args,cwd=self.root,capture_output=True,text=True,timeout=60)
        self.assertNotEqual(result.returncode,0,result.stderr)
        if message:
            self.assertIn(message,result.stderr)
        return result

    def launch(self, args):
        self.serial += 1
        log = self.root/('process-'+str(self.serial)+'.log')
        with log.open('w') as handle:
            process = subprocess.Popen([str(BINARY)]+args,cwd=self.root,stdout=subprocess.DEVNULL,
                stderr=handle,start_new_session=True)
        self.children.append(process)
        return process,log

    def wait_log(self, process, log, text, regex=False):
        deadline = time.monotonic()+25
        while time.monotonic() < deadline:
            content = log.read_text()
            if re.search(text,content) if regex else text in content:
                return content
            if process.poll() is not None:
                self.fail('Process ended before '+text+':\n'+content)
            time.sleep(0.005)
        self.fail('Timed out waiting for '+text+':\n'+log.read_text())

    def interrupt(self, process, log, signum=signal.SIGUSR1):
        if signum == signal.SIGKILL:
            os.killpg(process.pid,signum)
        else:
            process.send_signal(signum)
        code = process.wait(timeout=25)
        self.assertEqual(code,-signal.SIGKILL if signum==signal.SIGKILL else 75,log.read_text())

    def checkpoint_digest(self, output):
        return hashlib.sha256((output/'checkpoint/state.bin').read_bytes()).hexdigest()

    def manifest(self, output):
        with (output/'manifest.tsv').open() as handle:
            return {r['key']:r['value'] for r in csv.DictReader(handle,delimiter='\t')}

    def compare(self, expected, actual):
        self.assertTrue((actual/'SUCCESS').exists())
        for name in ('markers.fasta','markers.tsv','target_taxids.txt','target_taxonomy.tsv'):
            self.assertEqual((expected/name).read_bytes(),(actual/name).read_bytes(),name)
        a,b = self.manifest(expected),self.manifest(actual)
        for key in ('target_defline_rows','target_bases','candidate_kmers','all_defline_rows',
                    'background_defline_rows','background_bases','removed_kmers','retained_kmers'):
            self.assertEqual(a[key],b[key],key)
        self.assertEqual(b['status'],'complete')
        self.assertGreaterEqual(int(b['resume_count']),1)

    def paused_background(self, name='paused', signum=signal.SIGUSR1):
        output = self.root/name
        self.configure(row_delay=0.28)
        process,log = self.launch(self.args(output)+['-t','4'])
        self.wait_log(process,log,r'Saved background scan after [1-9][0-9]* records',regex=True)
        self.interrupt(process,log,signum)
        self.assertFalse((output/'SUCCESS').exists())
        return output

    def test_target_loading_sigterm_resume_all_k_boundaries(self):
        for k in (1,32,33,42,64):
            with self.subTest(k=k):
                self.configure()
                baseline=self.root/('baseline-'+str(k))
                self.run_ok(self.args(baseline,k))
                self.configure(target_row_delay=0.28)
                output=self.root/('target-'+str(k))
                process,log=self.launch(self.args(output,k))
                self.wait_log(process,log,r'Saved target loading after [1-9][0-9]* records',regex=True)
                self.interrupt(process,log,signal.SIGTERM)
                self.assertEqual(self.manifest(output)['status'],'paused')
                self.configure()
                result=self.run_ok(['-resume',str(output),'-t','8'])
                self.assertIn('Saved prefix verified',result.stderr)
                self.compare(baseline,output)

    def test_background_sigkill_resume_and_repeated_pause_during_replay(self):
        baseline=self.root/'baseline'
        self.run_ok(self.args(baseline))
        output=self.paused_background(signum=signal.SIGKILL)
        digest=self.checkpoint_digest(output)
        self.configure(row_delay=0.28)
        process,log=self.launch(['-resume',str(output)])
        self.wait_log(process,log,'Replay saved input prefix')
        self.interrupt(process,log,signal.SIGINT)
        self.assertEqual(self.checkpoint_digest(output),digest)
        self.configure()
        self.run_ok(['-resume',str(output),'-t','30'])
        self.compare(baseline,output)

    def test_saved_exclusions_across_chunk_boundaries(self):
        # Exclusion flags must survive checkpointing, with worker queues drained.
        prefix='N'*(1024*1024-20)
        self.background=''.join(row('b'+str(i),prefix+s+'N'*100,'9606') for i,s in enumerate(self.sequences[:6]))
        (self.db/'all.tsv').write_text(self.background)
        baseline=self.root/'baseline'
        self.run_ok(self.args(baseline)+['-t','4'])
        output=self.paused_background()
        self.configure()
        self.run_ok(['-resume',str(output),'-t','1'])
        self.compare(baseline,output)
        self.assertGreater(int(self.manifest(output)['removed_kmers']),0)

    @unittest.skipUnless(GZIP_BINARY,'gzip required')
    def test_gzip_fasta_target_and_background_resume_with_compressed_map_and_dumps(self):
        content=''.join('>s'+str(i)+'\n'+s+'\n' for i,s in enumerate(self.sequences))
        mapping=''.join('s'+str(i)+'\t'+('10226' if i<4 else '9606')+'\n' for i in range(8))
        fasta=self.root/'references.fa.gz'
        fasta.write_bytes(gzip.compress(content[:47].encode())+gzip.compress(content[47:].encode()))
        paths={'-fasta':fasta}
        for flag,name,data in [('-taxidmap','mapping.tsv',mapping.encode()),
                               ('-taxnodes','nodes.dmp',(self.root/'nodes.dmp').read_bytes()),
                               ('-taxnames','names.dmp',(self.root/'names.dmp').read_bytes())]:
            path=self.root/(name+'.gz'); path.write_bytes(gzip.compress(data)); paths[flag]=path
        reader=self.root/'slow-gzip.py'
        counter=self.root/'gzip-count'
        # Only the selected FASTA pass is slow. Mapping and taxonomy retain normal behavior.
        def set_reader(slow_pass):
            reader.write_text('#!'+sys.executable+'\n'+
                'import gzip,pathlib,sys,time\npath=pathlib.Path(sys.argv[4])\n'+
                'counter=pathlib.Path('+repr(str(counter))+')\n'+
                'n=0\nif path.name=="references.fa.gz":\n'+
                ' n=int(counter.read_text())+1 if counter.exists() else 1\n counter.write_text(str(n))\n'+
                'for line in gzip.decompress(path.read_bytes()).splitlines(keepends=True):\n'+
                ' sys.stdout.buffer.write(line); sys.stdout.buffer.flush()\n'+
                ' if n=='+str(slow_pass)+': time.sleep(0.16)\n')
            reader.chmod(0o755)
        common=['-resume-mode','replay','-taxid','10226','-checkpoint-seconds','1']
        for flag,path in paths.items(): common += [flag,str(path)]
        baseline=self.root/'baseline'
        self.run_ok(common+['-gzip',GZIP_BINARY,'-out',str(baseline)])
        for phase,slow_pass in [('target loading',1),('background scan',2)]:
            counter.write_text('0'); set_reader(slow_pass)
            output=self.root/('gzip-'+str(slow_pass))
            process,log=self.launch(common+['-gzip',str(reader),'-out',str(output)])
            self.wait_log(process,log,'Saved '+phase+r' after [1-9][0-9]* records',regex=True)
            self.interrupt(process,log)
            self.run_ok(['-resume',str(output),'-gzip',GZIP_BINARY,'-t','4'])
            self.compare(baseline,output)
            if slow_pass == 2:
                # An output checkpoint must still reject changed mapping content.
                (output/'SUCCESS').unlink()
                paths['-taxidmap'].write_bytes(gzip.compress(mapping.replace('10226','9606',1).encode()))
                self.run_bad(['-resume',str(output),'-gzip',GZIP_BINARY],'Checkpoint inputs changed')

    def test_changed_prefix_metadata_taxonomy_and_truncation_rejected(self):
        output=self.paused_background()
        digest=self.checkpoint_digest(output)
        original=(self.db/'all.tsv').read_text()
        nodes=(self.root/'nodes.dmp').read_bytes()
        cases=['prefix','truncated','metadata','taxonomy']
        for case in cases:
            with self.subTest(case=case):
                self.configure()
                (self.db/'all.tsv').write_text(original)
                (self.root/'nodes.dmp').write_bytes(nodes)
                if case=='prefix': (self.db/'all.tsv').write_text(original.replace('t0','altered',1))
                if case=='truncated': (self.db/'all.tsv').write_text('')
                if case=='metadata': self.configure(database_info='changed database snapshot')
                if case=='taxonomy': (self.root/'nodes.dmp').write_bytes(nodes+b'999\t|\t1\t|\n')
                self.run_bad(['-resume',str(output)])
                self.assertFalse((output/'SUCCESS').exists())
                self.assertEqual(self.checkpoint_digest(output),digest)

    def test_concurrent_resume_lock_and_completed_run_protection(self):
        output=self.paused_background()
        digest=self.checkpoint_digest(output)
        self.configure(row_delay=0.3)
        process,log=self.launch(['-resume',str(output)])
        self.wait_log(process,log,'Replay saved input prefix')
        self.run_bad(['-resume',str(output)],'locked')
        self.assertEqual(self.checkpoint_digest(output),digest)
        self.interrupt(process,log)
        self.configure()
        self.run_ok(['-resume',str(output)])
        old=(output/'markers.tsv').read_bytes()
        self.run_bad(['-resume',str(output)],'already complete')
        self.assertEqual((output/'markers.tsv').read_bytes(),old)

    def test_output_restart_corruption_atomic_write_failure_and_partial_file(self):
        baseline=self.root/'baseline'
        self.run_ok(self.args(baseline))
        for case in ('output','truncated','corrupt','base','missing_base','partial','write_failure'):
            with self.subTest(case=case):
                output=self.root/case
                shutil.copytree(baseline,output)
                (output/'SUCCESS').unlink() # Simulate termination after the output checkpoint.
                state=output/'checkpoint/state.bin'
                base=output/'checkpoint/candidates.bin'
                if case=='truncated': state.write_bytes(state.read_bytes()[:-4])
                if case=='corrupt':
                    data=bytearray(state.read_bytes()); data[-24]^=1; state.write_bytes(data)
                if case=='base':
                    data=bytearray(base.read_bytes()); data[-24]^=1; base.write_bytes(data)
                if case=='missing_base': base.unlink()
                if case=='partial': (output/'checkpoint/state.bin.partial').write_bytes(b'interrupted write')
                if case=='write_failure': (output/'checkpoint/state.bin.partial').mkdir()
                digest=self.checkpoint_digest(output)
                if case in ('truncated','corrupt','base','missing_base','write_failure'):
                    self.run_bad(['-resume',str(output)])
                    self.assertFalse((output/'SUCCESS').exists())
                    self.assertEqual(self.checkpoint_digest(output),digest)
                    if case!='write_failure': continue
                    (output/'checkpoint/state.bin.partial').rmdir()
                result=self.run_ok(['-resume',str(output)])
                self.assertIn('Background scan already complete',result.stderr)
                self.compare(baseline,output)

    def test_plain_fasta_output_resume_and_input_binding(self):
        fasta=self.root/'references.fa'
        content=(ROOT/'examples/references_taxids.fasta').read_bytes()
        fasta.write_bytes(content)
        baseline=self.root/'baseline'
        args=['-fasta',str(fasta),'-taxnodes',str(self.root/'nodes.dmp'),'-taxnames',str(self.root/'names.dmp'),
              '-taxid','10226','-out',str(baseline)]
        self.run_ok(args)
        output=self.root/'resumed'; shutil.copytree(baseline,output); (output/'SUCCESS').unlink()
        self.run_ok(['-resume',str(output)])
        self.compare(baseline,output)
        (output/'SUCCESS').unlink()
        fasta.write_bytes(content.replace(b'A',b'C',1))
        self.run_bad(['-resume',str(output)],'Checkpoint inputs changed')

    def test_checkpoint_cli_and_disabled_mode(self):
        for args in (['-resume',str(self.root/'missing')],['-checkpoint-seconds','-1'],
                     ['-checkpoint-seconds','86401'],['-checkpoint-seconds','1.5'],
                     ['-resume',str(self.root),'-k','42'],['-resume',str(self.root),'-out','elsewhere']):
            self.run_bad(args)
        output=self.root/'disabled'
        args=self.args(output); args[-1]='0'
        self.run_ok(args)
        self.assertFalse((output/'checkpoint').exists())

    def test_slurm_wrapper_forwards_signal_and_waits_for_save(self):
        fake=self.root/'fake-clademers.py'
        output=self.root/'result'
        ready=self.root/'ready'
        recorded=self.root/'arguments.json'
        fake.write_text('#!'+sys.executable+'\n'+
            'import json,pathlib,signal,sys,time\n'+
            'result=pathlib.Path('+repr(str(output))+')\n'+
            'pathlib.Path('+repr(str(recorded))+').write_text(json.dumps(sys.argv[1:]))\n'+
            'if "-resume" in sys.argv:\n (result/"SUCCESS").write_text("complete")\n sys.exit(0)\n'+
            'def pause(signum,frame):\n time.sleep(0.15)\n'+
            ' (result/"checkpoint").mkdir(parents=True)\n'+
            ' (result/"checkpoint/state.bin").write_text("saved")\n sys.exit(75)\n'+
            'signal.signal(signal.SIGUSR1,pause)\n'+
            'pathlib.Path('+repr(str(ready))+').write_text("ready")\n'+
            'while True: time.sleep(1)\n')
        fake.chmod(0o755)
        script=self.root/'job.sh'
        text=(ROOT/'examples/nt_checkpoint.slurm').read_text()
        text=text.replace('CLADEMERS=/path/to/clademers','CLADEMERS='+str(fake))
        text=text.replace('RESULT=/scratch/yourname/placozoa_k42','RESULT='+str(output))
        script.write_text(text)
        subprocess.run(['bash','-n',str(script)],check=True)
        process=subprocess.Popen(['bash',str(script)],cwd=self.root,stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,text=True,start_new_session=True)
        self.children.append(process)
        deadline=time.monotonic()+10
        while not ready.exists() and process.poll() is None and time.monotonic()<deadline:
            time.sleep(0.005)
        self.assertTrue(ready.exists())
        process.send_signal(signal.SIGUSR1)
        stdout,stderr=process.communicate(timeout=10)
        self.assertEqual(process.returncode,75,stderr)
        self.assertTrue((output/'checkpoint/state.bin').exists())
        self.assertIn('resubmit',stdout)
        result=subprocess.run(['bash',str(script)],cwd=self.root,capture_output=True,text=True,timeout=10)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertIn('-resume',json.loads(recorded.read_text()))
        self.assertTrue((output/'SUCCESS').exists())

    @unittest.skipUnless(os.environ.get('BLAST_BIN') and SQLITE3_BINARY,'BLAST+ and SQLite required')
    def test_real_blast_and_sqlite_taxonomy_resume(self):
        bindir=Path(os.environ['BLAST_BIN'])
        real=bindir/'blastdbcmd'
        fasta=self.root/'references.fasta'; taxmap=self.root/'taxids.tsv'
        entries=[('b'+str(i),s,9606) for i,s in enumerate(self.sequences[:4])]
        entries += [('t'+str(i),s,10226) for i,s in enumerate(self.sequences)]
        fasta.write_text(''.join('>'+name+'\n'+seq+'\n' for name,seq,_ in entries))
        taxmap.write_text(''.join(name+'\t'+str(tax)+'\n' for name,_,tax in entries))
        database=self.root/'real-db'
        made=subprocess.run([str(bindir/'makeblastdb'),'-in',str(fasta),'-dbtype','nucl',
            '-parse_seqids','-taxid_map',str(taxmap),'-out',str(database)],capture_output=True,text=True)
        self.assertEqual(made.returncode,0,made.stderr)
        taxonomy=write_blast_taxonomy(self.root/'blast-taxonomy',self.root)
        common=['-resume-mode','replay','-blastdb',str(database),'-blasttax',str(taxonomy),'-sqlite3',SQLITE3_BINARY,
                '-taxid','10226','-checkpoint-seconds','1','-t','4']
        baseline=self.root/'baseline'
        self.run_ok(common+['-blastdbcmd',str(real),'-out',str(baseline)])
        wrapper=self.root/'slow-blast.py'
        wrapper.write_text('#!'+sys.executable+'\n'+
            'import subprocess,sys,time\n'+
            'p=subprocess.Popen(['+repr(str(real))+']+sys.argv[1:],stdout=subprocess.PIPE)\n'+
            'for line in p.stdout:\n sys.stdout.buffer.write(line); sys.stdout.buffer.flush()\n'+
            ' if "-entry" in sys.argv: time.sleep(0.28)\n'+
            'sys.exit(p.wait())\n')
        wrapper.chmod(0o755)
        output=self.root/'resumed'
        process,log=self.launch(common+['-blastdbcmd',str(wrapper),'-out',str(output)])
        self.wait_log(process,log,r'Saved background scan after [1-9][0-9]* records',regex=True)
        self.interrupt(process,log)
        self.run_ok(['-resume',str(output),'-blastdbcmd',str(real),'-t','1'])
        self.compare(baseline,output)

    @unittest.skipUnless(os.environ.get('CLADEMERS_RESUME_BINARY'),'optional second-architecture binary')
    def test_resume_on_second_architecture(self):
        other=Path(os.environ['CLADEMERS_RESUME_BINARY'])
        baseline=self.root/'baseline'
        self.run_ok(self.args(baseline,k=64))
        self.configure(target_row_delay=0.28)
        output=self.root/'target'
        process,log=self.launch(self.args(output,k=64))
        self.wait_log(process,log,r'Saved target loading after [1-9][0-9]* records',regex=True)
        self.interrupt(process,log)
        self.configure()
        self.run_ok(['-resume',str(output),'-t','4'],binary=other)
        self.compare(baseline,output)
        baseline2=self.root/'baseline2'
        self.run_ok(self.args(baseline2))
        output=self.paused_background(name='background')
        self.configure()
        self.run_ok(['-resume',str(output),'-t','8'],binary=other)
        self.compare(baseline2,output)


if __name__ == '__main__':
    unittest.main(verbosity=2)
