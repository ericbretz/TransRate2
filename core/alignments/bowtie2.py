import os
import subprocess
import sys
import time
from core.utils.logging import LoggingSubprocess

class Bowtie2:
    def __init__(self, dict_dir, dict_file, dict_info, printout_class, transrate_logger=None):
        self.printout_class     = printout_class
        self.printout           = printout_class.printout
        self.transrate_logger   = transrate_logger
        self.logging_subprocess = LoggingSubprocess(transrate_logger, 'bowtie2') if transrate_logger else None
        self.dict_dir           = dict_dir
        self.dict_file          = dict_file
        self.dict_info          = dict_info
        self.threads            = self.dict_info['threads']
        self.mode               = self.dict_info['mode']
        self.assembly_path      = self.dict_file['assembly']
        self.aligner_dir        = self.dict_dir['temp_aligner']
        self.aligner_index      = self.dict_dir['temp_aligner_index']
        self.aligner_prefix     = self.dict_file['aligner_prefix']

        self.printout('subtitle', 'Bowtie2')

        return
    
    def index(self):
        self.printout_class.start_progress_stage("Indexing")
        
        self.printout_class.update_progress_bar(1, 3, "Preparing")
        time.sleep(0.1)
        
        self.printout_class.update_progress_bar(2, 3, "Building index", {"input_file": os.path.basename(str(self.assembly_path))})
        
        bt2_index_cmd = ['bowtie2-build', '--threads', str(self.threads), str(self.assembly_path), str(self.aligner_prefix)]
        
        if self.logging_subprocess:
            returncode, stdout, stderr = self.logging_subprocess.run_with_logging(bt2_index_cmd, "indexing")
        else:
            bt2_index_run  = subprocess.Popen(bt2_index_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
            stdout, stderr = bt2_index_run.communicate()
            returncode     = bt2_index_run.returncode
            stdout         = stdout.decode('utf-8') if stdout else ""
            stderr         = stderr.decode('utf-8') if stderr else ""
        
        if returncode != 0:
            self.printout('error', f"Bowtie2 indexing failed: {stderr}")
            sys.exit(1)
        
        index_size = self._get_index_size()
        self.printout_class.update_progress_bar(3, 3, "Finalizing", {"Index Size": index_size})
        
        final_metrics = {
            "# Index Files": len([f for f in os.listdir(self.aligner_index) if f.startswith(os.path.basename(str(self.aligner_prefix)))]),
            "Index Size": index_size
        }
        
        print()
        for key, value in final_metrics.items():
            print(f"  {key:<25} {value}")
        print()
        
        return returncode
    
    def align(self):
        self.printout_class.start_progress_stage("Alignment")
        
        mode_type = "single-end" if self.dict_info['mode'] == 1 else "paired-end"
        self.printout_class.update_progress_bar(1, 3, "Preparing")
        time.sleep(0.1)
        
        self.printout_class.update_progress_bar(2, 3, "Aligning reads", {"sensitivity": "very-sensitive-local"})
        
        # bt2_cmd = ['bowtie2', '--threads', str(self.threads), '--very-sensitive', '--phred33', '--no-unal', '--rdg', '1000,1000', '--rfg', '1000,1000', '-x', str(self.aligner_prefix)]
        bt2_cmd = ['bowtie2', '--threads', str(self.threads), '--very-sensitive', '-x', str(self.aligner_prefix)]
        bt2_cmd.extend(['-U', self.dict_file['single'], '-S', self.dict_file['aligner_bam']] if self.dict_info['mode'] == 1 else ['-1', self.dict_file['left'], '-2', self.dict_file['right'], '-S', self.dict_file['aligner_bam']])
        
        if self.logging_subprocess:
            returncode, stdout, stderr = self.logging_subprocess.run_with_logging(bt2_cmd, "alignment")
        else:
            bt2_run        = subprocess.Popen(bt2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
            stdout, stderr = bt2_run.communicate()
            returncode     = bt2_run.returncode
            stdout         = stdout.decode('utf-8') if stdout else ""
            stderr         = stderr.decode('utf-8') if stderr else ""
        
        if returncode != 0:
            self.printout('error', f"Bowtie2 alignment failed: {stderr}")
            sys.exit(1)
        
        output_size     = self._get_bam_size()
        alignment_stats = self._parse_alignment_stats(stderr)
        self.printout_class.update_progress_bar(3, 3, "Processing output", {"Size": output_size})
        
        final_metrics = {
            "Size": output_size,
        }
        final_metrics.update(alignment_stats)
        
        print()
        for key, value in final_metrics.items():
            print(f"  {key:<25} {value}")
        print()
        
        return returncode
    
    def _get_index_size(self):
        try:
            total_size = 0
            for f in os.listdir(self.aligner_index):
                if f.startswith(os.path.basename(str(self.aligner_prefix))):
                    file_path = os.path.join(self.aligner_index, f)
                    total_size += os.path.getsize(file_path)
            return f"{total_size / (1024*1024):.1f}MB"
        except:
            return "unknown"
    
    def _get_bam_size(self):
        try:
            size = os.path.getsize(self.dict_file['aligner_bam'])
            return f"{size / (1024*1024):.1f}MB"
        except:
            return "unknown"
    
    def _parse_alignment_stats(self, stderr_output):
        stats = {}
        try:
            lines = stderr_output.split('\n')
            for line in lines:
                if 'reads; of these:' in line:
                    stats['Total Reads'] = line.split()[0]
                elif 'overall alignment rate' in line:
                    stats['Alignment Rate'] = line.split()[0]
        except:
            pass
        return stats