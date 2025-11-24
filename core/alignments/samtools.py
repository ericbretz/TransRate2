import os
import subprocess
import sys
import time
from core.utils.logging import LoggingSubprocess

class Samtools:
    def __init__(self, dict_dir, dict_file, dict_info, printout_class, transrate_logger=None):
        self.printout_class     = printout_class
        self.printout           = printout_class.printout
        self.transrate_logger   = transrate_logger
        self.logging_subprocess = LoggingSubprocess(transrate_logger, 'samtools') if transrate_logger else None
        self.dict_dir           = dict_dir
        self.dict_file          = dict_file
        self.dict_info          = dict_info
        self.threads            = self.dict_info['threads']
        self.mode               = self.dict_info['mode']
        
        self.printout('subtitle', 'Samtools')

        return
    
    def sort(self):
        self.printout_class.start_progress_stage("Sorting")
        
        input_size = self._get_input_size()
        self.printout_class.update_progress_bar(1, 2, "Preparing")
        time.sleep(0.1)
        
        self.printout_class.update_progress_bar(2, 2, "Sorting BAM", {"input": os.path.basename(str(self.dict_file['salmon_bam']))})
        
        samtools_sort_cmd = ['samtools', 'sort', f'-@{self.threads}', self.dict_file['salmon_bam'], '-o', self.dict_file['samtools_bam']]
        
        if self.logging_subprocess:
            returncode, stdout, stderr = self.logging_subprocess.run_with_logging(samtools_sort_cmd, "sort")
        else:
            samtools_sort_run = subprocess.Popen(samtools_sort_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
            stdout, stderr    = samtools_sort_run.communicate()
            returncode        = samtools_sort_run.returncode
            stdout            = stdout.decode('utf-8') if stdout else ""
            stderr            = stderr.decode('utf-8') if stderr else ""
        
        if returncode != 0:
            self.printout('error', f"Samtools sort failed: {stderr}")
            sys.exit(1)
        
        output_size = self._get_output_size()
        
        final_metrics = {
            "Input Size" : input_size,
            "Output Size": output_size,
        }
        
        print()
        for key, value in final_metrics.items():
            print(f"  {key:<25} {value}")
        print()
        
        return returncode

    def index(self):
        self.printout_class.start_progress_stage("Indexing")
        
        bam_size = self._get_bam_size()
        self.printout_class.update_progress_bar(1, 2, "Preparing")
        time.sleep(0.1)
        
        self.printout_class.update_progress_bar(2, 2, "Creating index", {"input": os.path.basename(str(self.dict_file['samtools_bam']))})
        
        samtools_index_cmd = ['samtools', 'index', '-b', f'-@{self.threads}', str(self.dict_file['samtools_bam'])]
        
        if self.logging_subprocess:
            returncode, stdout, stderr = self.logging_subprocess.run_with_logging(samtools_index_cmd, "index")
        else:
            samtools_index_run = subprocess.Popen(samtools_index_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
            stdout, stderr     = samtools_index_run.communicate()
            returncode         = samtools_index_run.returncode
            stdout             = stdout.decode('utf-8') if stdout else ""
            stderr             = stderr.decode('utf-8') if stderr else ""
        
        if returncode != 0:
            self.printout('error', f"Samtools index failed: {stderr}")
            sys.exit(1)
        
        index_size = self._get_index_size()
        
        final_metrics = {
            "Bam Size"  : bam_size,
            "Index Size": index_size,
        }
        
        print()
        for key, value in final_metrics.items():
            print(f"  {key:<25} {value}")
        print()
        
        return returncode
    
    def _get_input_size(self):
        try:
            size = os.path.getsize(self.dict_file['salmon_bam'])
            return f"{size / (1024*1024):.1f}MB"
        except:
            return "unknown"
    
    def _get_output_size(self):
        try:
            size = os.path.getsize(self.dict_file['samtools_bam'])
            return f"{size / (1024*1024):.1f}MB"
        except:
            return "unknown"
    
    def _get_bam_size(self):
        try:
            size = os.path.getsize(self.dict_file['samtools_bam'])
            return f"{size / (1024*1024):.1f}MB"
        except:
            return "unknown"
    
    def _get_index_size(self):
        try:
            index_file = str(self.dict_file['samtools_bam']) + '.bai'
            if os.path.exists(index_file):
                size = os.path.getsize(index_file)
                return f"{size / (1024*1024):.1f}MB"
            return "created"
        except:
            return "unknown"
    
