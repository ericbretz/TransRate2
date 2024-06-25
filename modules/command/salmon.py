import os
import subprocess
import sys

class Salmon:

    def __init__(self, main):
        pass

    def salmon(self, main):
        # Logging Entry
        main.LOG.log_time(main, 'salmon', 'start')
        main.LOG.log_set(main, 'salmon')
        main.STAGE = 'Salmon'

        salmon_cmd      = ['salmon', 'quant', '--libType', 'IU', '--alignments', main.BAM_ALIGNER, '--targets', main.ASSEMBLY, f'--threads={main.THREADS}', '--sampleOut', '--sampleUnaligned', '--output', main.SALMON_PATH]
        salmon_run = subprocess.Popen(salmon_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr = salmon_run.communicate()
        returncode = salmon_run.returncode

        # Logging Exit
        main.LOG.log_time(main, 'salmon', 'end')
        main.LOG.log_write(main, stdout, stderr)

        # Error Handling
        if salmon_run.returncode != 0:
            main.ERROR = True
            # while not main.STOPPED:
            #     pass
            # print('\n')
            # print('Salmon Failed')
            # print(f'Check Logs at {main.LOG_PATH}')
            # try:
            #     main.OUTPUTTHREAD.join()
            # except:
            #     pass
            # # sys.exit(1)
            main.LOG.error_out(main, 'Salmon', 'Salmon failed')

        if os.path.exists(os.path.join(main.SALMON_PATH, 'postSample.bam')):
            new_name = os.path.join(main.SALMON_PATH, f'{main.ASSEMBLY_NAME}.postSample.bam')
            os.rename(os.path.join(main.SALMON_PATH, 'postSample.bam'), new_name)
            main.BAM_SALMON = new_name
        else:
            raise Exception('Error: Salmon failed to generate postSample.bam file')
        


        # if os.path.exists(os.path.join(main.SALMON_PATH, 'quant.sf')):
        #     main.SALMON_QUANT = os.path.join(main.SALMON_PATH, 'quant.sf')
        # else:
        #     raise Exception('Error: Salmon failed to generate quant.sf file')

