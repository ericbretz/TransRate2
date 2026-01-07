import os
import subprocess
import sys
import time
import json
from core.utils.logging import LoggingSubprocess

class Salmon:
    def __init__(self, dict_dir, dict_file, dict_info, printout_class, transrate_logger=None):
        self.printout_class     = printout_class
        self.printout           = printout_class.printout
        self.transrate_logger   = transrate_logger
        self.logging_subprocess = LoggingSubprocess(transrate_logger, 'salmon') if transrate_logger else None
        self.dict_dir           = dict_dir
        self.dict_file          = dict_file
        self.dict_info          = dict_info
        self.threads            = self.dict_info['threads']
        self.mode               = self.dict_info['mode']
        self.printout('subtitle', 'Salmon')
        return
    
    def quant(self):
        self.printout_class.start_progress_stage("Quantification")
        
        self.printout_class.update_progress_bar(1, 3, "Preparing")
        time.sleep(0.1)
        
        self.printout_class.update_progress_bar(2, 3, "Quantifying", {"input": os.path.basename(str(self.dict_file['aligner_bam']))})
        
        salmon_cmd = ['salmon', 'quant', '--libType', 'A', '--alignments', str(self.dict_file['aligner_bam']), '--targets', str(self.dict_file['assembly']), f'--threads={self.threads}', '--sampleOut', '--sampleUnaligned', '--output', str(self.dict_dir['temp_salmon'])]
        
        if self.logging_subprocess:
            returncode, stdout, stderr = self.logging_subprocess.run_with_logging(salmon_cmd, "quantification")
        else:
            salmon_run     = subprocess.Popen(salmon_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
            stdout, stderr = salmon_run.communicate()
            returncode     = salmon_run.returncode
            stdout         = stdout.decode('utf-8') if stdout else ""
            stderr         = stderr.decode('utf-8') if stderr else ""
        
        if returncode != 0:
            self.printout('error', f"Salmon quantification failed: {stderr}")
            sys.exit(1)
        
        output_size = self._get_output_size()
        log_stats   = self._parse_salmon_log()
        meta_stats  = self._parse_meta_info()
        quant_stats = self._parse_quant_file()
        self.printout_class.update_progress_bar(3, 3, "Processing output", {"Size": output_size})
        
        if os.path.exists(os.path.join(self.dict_dir['temp_salmon'], 'postSample.bam')):
            os.rename(os.path.join(self.dict_dir['temp_salmon'], 'postSample.bam'), self.dict_file['salmon_bam'])
        
        final_metrics = {
            "Size": output_size,
        }
        final_metrics.update(log_stats)
        final_metrics.update(meta_stats)
        final_metrics.update(quant_stats)
        
        print()
        for key, value in final_metrics.items():
            print(f"  {key:<25} {value}")
        print()
        
        return returncode
    
    def _get_output_size(self):
        try:
            if os.path.exists(self.dict_file['salmon_bam']):
                size = os.path.getsize(self.dict_file['salmon_bam'])
                return f"{size / (1024*1024):.1f}MB"
            temp_bam = os.path.join(self.dict_dir['temp_salmon'], 'postSample.bam')
            if os.path.exists(temp_bam):
                size = os.path.getsize(temp_bam)
                return f"{size / (1024*1024):.1f}MB"
            return "unknown"
        except:
            return "unknown"
    
    def _parse_salmon_log(self):
        stats = {}
        log_file = os.path.join(self.dict_dir['temp_salmon'], 'logs', 'salmon_quant.log')
        try:
            if os.path.exists(log_file):
                with open(log_file, 'r') as f:
                    content = f.read()
                    lines = content.split('\n')
                    for line in lines:
                        if 'Total # of mapped reads :' in line:
                            stats['Mapped Reads'] = line.split(':')[1].strip()
                        elif '# of uniquely mapped reads :' in line:
                            stats['Unique Reads'] = line.split(':')[1].strip()
                        elif '# ambiguously mapped reads :' in line:
                            stats['Ambiguous Reads'] = line.split(':')[1].strip()
                        elif 'Computed' in line and 'rich equivalence classes' in line:
                            parts = line.split()
                            for i, part in enumerate(parts):
                                if part.replace(',', '').isdigit():
                                    stats['Equiv Classes'] = part.replace(',', '')
                                    break
        except:
            pass
        return stats
    
    def _parse_meta_info(self):
        stats = {}
        meta_file = os.path.join(self.dict_dir['temp_salmon'], 'aux_info', 'meta_info.json')
        try:
            if os.path.exists(meta_file):
                with open(meta_file, 'r') as f:
                    meta = json.load(f)
                    stats['Processed Fragments'] = f"{meta.get('num_processed', 0):,}"
                    stats['Mapping Rate']        = f"{meta.get('percent_mapped', 0):.2f}%"
                    stats['Valid Targets']       = f"{meta.get('num_valid_targets', 0):,}"
                    stats['Fragment Length']     = f"{meta.get('frag_length_mean', 0):.1f}"
                    stats['Library Type']        = meta.get('library_types', ['unknown'])[0] if meta.get('library_types') else 'unknown'
        except:
            pass
        return stats
    
    def _parse_quant_file(self):
        stats = {}
        quant_file = os.path.join(self.dict_dir['temp_salmon'], 'quant.sf')
        try:
            if os.path.exists(quant_file):
                with open(quant_file, 'r') as f:
                    lines                = f.readlines()
                    transcript_count     = len(lines) - 1
                    stats['Transcripts'] = f"{transcript_count:,}"
                    
                    total_tpm = 0
                    expressed_transcripts = 0
                    for line in lines[1:]:
                        parts = line.strip().split('\t')
                        if len(parts) >= 4:
                            tpm = float(parts[3])
                            total_tpm += tpm
                            if tpm > 1.0:
                                expressed_transcripts += 1
                    
                    stats['Expressed (TPM>1)'] = f"{expressed_transcripts:,}"
                    stats['Expression Coverage'] = f"{(expressed_transcripts/transcript_count*100):.1f}%" if transcript_count > 0 else "0%"
        except:
            pass
        return stats