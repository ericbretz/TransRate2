import os
import subprocess
import sys

class SNAP:
    def __init__(self, main):
        pass

    def snap_index(self, main):
        # Logging Entry
        main.LOG.log_time(main, 'snap_index', 'start')
        main.LOG.log_set(main, 'snap_index')
        main.STAGE = 'Snap Index'

        # Create snap index
        snap_index_cmd = ['snap-aligner', 'index', main.ASSEMBLY, main.ALIGNER_INDEX, '-s', '23', f'-t{main.THREADS}', '-bSpace', '-locationSize', '4']
        snap_index_run = subprocess.Popen(snap_index_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)

        stdout, stderr = snap_index_run.communicate()
        returncode = snap_index_run.returncode

        # Logging Exit
        main.LOG.log_time(main, 'snap_index', 'end')
        main.LOG.log_write(main, stdout, stderr)

        # Error Handling
        if snap_index_run.returncode != 0:
            main.ERROR = True

            try:
                main.OUT_THREAD.join()
            except:
                pass

            print('\n')
            print('Snap Index Failed')
            print(f'Check Logs at {main.LOG_PATH}')

            sys.exit(1)

    def snap(self, main):
        # Logging Entry
        main.LOG.log_time(main, 'snap_paired', 'start')
        main.LOG.log_set(main, 'snap_paired')
        main.STAGE = 'Snap'

        if main.READMODE == 2:
            snap_cmd    = ['snap-aligner', 'paired', main.ALIGNER_INDEX, main.LEFT, main.RIGHT, '-o', main.BAM_ALIGNER, '-s', '0', '1000', '-H', '300000',  '-h', '2000', '-d', '30', '-t', f'{main.THREADS}', '-b', '-M', '-D', '5', '-om', '5', '-omax', '10']
        else:
            snap_cmd    = ['snap-aligner', 'single', main.ALIGNER_INDEX, main.SINGLE, '-o', main.BAM_ALIGNER, '-h', '2000', '-d', '30', '-t', f'{main.THREADS}', '-b', '-D', '5', '-om', '5', '-omax', '10']
        snap_run        = subprocess.Popen(snap_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)

        stdout, stderr  = snap_run.communicate()
        returncode      = snap_run.returncode

        # Logging Exit
        main.LOG.log_time(main, 'snap_paired', 'end')
        main.LOG.log_write(main, stdout, stderr)

        # Get read count
        try:
            count = stdout.decode('utf-8').split('\n')
            if main.READMODE == 1:
                main.READ_COUNT = int(count[3].split(' ')[0].replace(',', ''))
            else:
                main.READ_COUNT = int(int(count[3].split(' ')[0].replace(',', '')) / 2)
        except:
            raise Exception('Error: Snap failed to align reads')
        
        # Error Handling
        if snap_run.returncode != 0:
            main.ERROR = True

            try:
                main.OUT_THREAD.join()
            except:
                pass

            print('\n')
            print('Snap Failed')
            print(f'Check Logs at {main.LOG_PATH}')

            sys.exit(1)
