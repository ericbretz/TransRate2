import os
import subprocess
import sys

class Samtools:
    def __init__(self, main):
        pass

    def samtools_sort(self, main):
        # Logging Entry
        main.LOG.log_time(main, 'samtools_sort', 'start')
        main.LOG.log_set(main, 'samtools_sort')
        main.STAGE = 'Samtools Sort'

        # Create sorted bam
        samtools_sort_cmd = ['samtools', 'sort', f'-@{main.THREADS}', main.BAM_SALMON, '-o', main.BAM_SORTED]
        samtools_sort_run = subprocess.Popen(samtools_sort_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr = samtools_sort_run.communicate()
        returncode = samtools_sort_run.returncode

        # Logging Exit
        main.LOG.log_time(main, 'samtools_sort', 'end')
        main.LOG.log_write(main, stdout, stderr)

        # Error Handling
        if samtools_sort_run.returncode != 0:
            main.ERROR = True
            # while not main.STOPPED:
            #     pass
            # print('\n')
            # print('Samtools Sort Failed')
            # print(f'Check Logs at {main.LOG_PATH}')
            # try:
            #     main.OUTPUTTHREAD.join()
            # except:
            #     pass
            # sys.exit(1)
            main.LOG.error_out(main, 'Samtools Sort', 'Samtools Sort Failed')

    def samtools_index(self, main):
        # Logging Entry
        main.LOG.log_time(main, 'samtools_index', 'start')
        main.LOG.log_set(main, 'samtools_index')
        main.STAGE = 'Samtools Index'

        samtools_index_cmd = ['samtools', 'index', '-b', f'-@{main.THREADS}', main.BAM_SORTED]
        samtools_index_run = subprocess.Popen(samtools_index_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr = samtools_index_run.communicate()

        # Logging Exit
        main.LOG.log_time(main, 'samtools_index', 'end')
        main.LOG.log_write(main, stdout, stderr)

        # Error Handling
        if samtools_index_run.returncode != 0:
            main.ERROR = True
            # while not main.STOPPED:
            #     pass
            # print('\n')
            # print('Samtools Index Failed')
            # print(f'Check Logs at {main.LOG_PATH}')
            # try:
            #     main.OUTPUTTHREAD.join()
            # except:
            #     pass
            # sys.exit(1)
            main.LOG.error_out(main, 'Samtools Index', 'Samtools Index Failed')
    
        # if os.path.exists(os.path.join(main.BAM_SALMON.strip('.bam') + '.sorted.bam')):
        #     main.BAM_SORTED = main.BAM_SALMON + '.sorted.bam'
        # else:
        #     print(os.path.join(main.BAM_SALMON.strip('.bam') + '.sorted.bam'))
        #     raise Exception('Error: Samtools failed to generate sorted bam file')
                    
        # if os.path.exists(main.BAM_SALMON + '.sorted.bam.bai'):
        #     main.BAM_INDEX = main.BAM_SALMON + '.sorted.bam.bai'
        # else:
        #     raise Exception('Error: Samtools failed to generate index file')