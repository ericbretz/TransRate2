import os
import subprocess
import sys

class BT2:
    def __init__(self, main):
        pass

    def bowtie2_index(self, main):
        # Logging Entry
        main.LOG.log_time(main, 'bowtie2_index', 'start')
        main.LOG.log_set(main, 'bowtie2_index')
        main.STAGE = 'Bowtie2 Index'

        # Create bowtie2 index
        bt2_index_cmd   = ['bowtie2-build', '--threads', str(main.THREADS), main.ASSEMBLY, main.ALIGNER_INDEX]
        bt2_index_run   = subprocess.Popen(bt2_index_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr  = bt2_index_run.communicate()
        returncode      = bt2_index_run.returncode

        # Logging Exit
        main.LOG.log_time(main, 'bowtie2_index', 'end')
        main.LOG.log_write(main, stdout, stderr)

        # Error Handling
        if bt2_index_run.returncode != 0:
            self.ERROR = True

            try:
                main.OUT_THREAD.join()
            except:
                pass

            print('\n')
            print('Bowtie2 Index Failed')
            print(f'Check Logs at {self.LOGDIR}')

            # sys.exit(1)
            main.LOG.error_out(main, 'Bowtie2 Index', 'Bowtie2 index failed')

    def bowtie2(self, main):
        # Logging Entry
        main.LOG.log_time(main, 'bowtie2', 'start')
        main.LOG.log_set(main, 'bowtie2')
        main.STAGE = 'Bowtie2'

        if main.READMODE == 2:
            bt2_cmd = ['bowtie2', '--threads', f'{main.THREADS}', '--very-sensitive', '--phred33', '--no-unal', '--no-mixed', '--no-discordant', '--rdg', '1000,1000', '--rfg', '1000,1000', '-x', main.ALIGNER_INDEX, '-1', main.LEFT, '-2', main.RIGHT, '-S', main.BAM_ALIGNER]
        else:
            bt2_cmd = ['bowtie2', '--threads', f'{main.THREADS}', '--very-sensitive', '--phred33', '--no-unal', '--no-mixed', '--no-discordant', '--rdg', '1000,1000', '--rfg', '1000,1000', '-x', main.ALIGNER_INDEX, '-U', main.SINGLE, '-S', main.BAM_ALIGNER]

        bt2_run = subprocess.Popen(bt2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr = bt2_run.communicate()
        returncode = bt2_run.returncode

        count = stderr.decode('utf-8').split(' ')[0]
        main.READ_COUNT = int(count)


        # Logging Exit
        main.LOG.log_time(main, 'bowtie2', 'end')
        main.LOG.log_write(main, stdout, stderr)

        # Error Handling
        if not main.READ_COUNT:
            # print('Error: Could not get read count from Bowtie2 output')
            # sys.exit(1)
            main.LOG.error_out(main, 'Bowtie2', 'Could not get read count from Bowtie2 output')

        if bt2_run.returncode != 0:
            main.ERROR = True

            try:
                main.OUT_THREAD.join()
            except:
                pass

            print('\n')
            print('Bowtie2 Failed')
            print(f'Check Logs at {main.LOG_PATH}')

            # sys.exit(1)
            main.LOG.error_out(main, 'Bowtie2', 'Bowtie2 failed')

        if os.path.exists(os.path.join(main.ALIGNER_PATH, main.ASSEMBLY_NAME + '.bam')):
            new_name = os.path.join(main.ALIGNER_PATH, main.ASSEMBLY_NAME + '.aligned.bam')
            os.rename(os.path.join(main.ALIGNER_PATH, main.ASSEMBLY_NAME + '.bam'), new_name)
            main.BAM_ALIGNER = new_name
        else:
            raise Exception('Error: Bowtie2 failed to create bam file')

