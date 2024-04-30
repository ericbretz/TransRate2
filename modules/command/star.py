import os
import subprocess
import sys
import psutil

class STAR:
    def __init__(self, main):
        pass

    def star_index(self, main):
        # Logging Entry
        main.LOG.log_time(main, 'star_index', 'start')
        main.LOG.log_set(main, 'star_index')
        main.STAGE = 'STAR Index'

        ram = psutil.virtual_memory().available

        # Create STAR index
        star_index_cmd = ['STAR', '--runThreadN', str(main.THREADS), '--runMode', 'genomeGenerate', '--genomeDir', main.ALIGNER_INDEX, '--genomeFastaFiles', main.ASSEMBLY, '--genomeSAindexNbases', '11', '--limitGenomeGenerateRAM', str(ram), '--outFileNamePrefix', os.path.join(main.ALIGNER_INDEX, 'star_')]
        star_index_run = subprocess.Popen(star_index_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr = star_index_run.communicate()
        returncode = star_index_run.returncode

        # Logging Exit
        main.LOG.log_time(main, 'star_index', 'end')
        main.LOG.log_write(main, stdout, stderr)

        if stderr:
            err = stderr.decode('utf-8')
            if '--limitGenomeGenerateRAM' in err:
                err = err.split('SOLUTION: please specify --limitGenomeGenerateRAM not less than ')[-1]
                err = err.split(' ')[0]
                err = f'{int(err) / (1024.0 ** 3):.2f} GB'
                ram = f'{ram / (1024.0 ** 3):.2f} GB'
                top   = f'\033[0;31m  ┌{"─" * 29}\033[m   \033[mSTAR Error   \033[0;31m{"─" * 29}┐\033[m'
                middle = f'Error: Need at least {err} RAM but only have {ram} RAM'
                middle = f'\033[0;31m  │{middle:^74}│\033[m'
                bottom= f'\033[0;31m  └{"─" * 74}┘\033[m'
                print(f'{top}\n{middle}\n{bottom}')
                main.ERROR = True
                sys.exit(1)

        # Error Handling
        if star_index_run.returncode != 0:
            self.ERROR = True

            try:
                main.OUT_THREAD.join()
            except:
                pass

            print('\n')
            print('STAR Index Failed')
            print(f'Check Logs at {self.LOGDIR}')

            sys.exit(1)

    def star(self, main):
        # Logging Entry
        main.LOG.log_time(main, 'star', 'start')
        main.LOG.log_set(main, 'star')
        main.STAGE = 'STAR'

        prefix = os.path.join(main.ALIGNER_PATH, main.ASSEMBLY_NAME + '_')

        reads = [main.SINGLE] if main.READMODE == 1 else [main.LEFT, main.RIGHT]
        star_cmd = ['STAR', '--runThreadN', str(main.THREADS), '--genomeDir', main.ALIGNER_INDEX, '--readFilesIn']
        star_cmd.extend(reads)
        star_cmd.extend(['--outFileNamePrefix', prefix])
        
        if reads[0].endswith('.gz'):
            star_cmd.extend(['--readFilesCommand', 'gunzip', '-c', '--outSAMtype BAM Unsorted'])
        else:
            star_cmd.extend(['--outSAMtype BAM Unsorted'])

        star_run = subprocess.Popen(star_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr = star_run.communicate()
        returncode = star_run.returncode

        # Logging Exit
        main.LOG.log_time(main, 'star', 'end')
        main.LOG.log_write(main, stdout, stderr)

        # Error Handling
        if star_run.returncode != 0:
            main.ERROR = True

            try:
                main.OUT_THREAD.join()
            except:
                pass

            print('\n')
            print('STAR Failed')
            print(f'Check Logs at {main.LOG_PATH}')
            print(stderr.decode('utf-8'))

            sys.exit(1)

        count_file = os.path.join(main.ALIGNER_PATH, main.ASSEMBLY_NAME + '_Log.final.out')
        with open(count_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if 'Number of input reads' in line:
                    main.READ_COUNT = int(line.split('\t')[-1])

            if not main.READ_COUNT:
                print('Error: Could not get read count from STAR output')
                sys.exit(1)