import time
import sys

class Output:
    def __init__(self, main):
        # Stage Info
        self.STAGE = None
        self.DESC = {
                'file'    : 'Preparing Database',
                'frag'    : 'Scanning for Fragmentation',
                'seq'     : 'Assessing Sequence Quality',
                'good'    : 'Determining Good/Bad Reads',
                'base'    : 'Locating Uncovered Bases',
                'sgmt'    : 'Calculating Segmentation',
                'plots'   : 'Creating Plots of Results',
                'csv'     : 'Writing Results to CSV',
                }
        
        # Spinner
        self.SPINCOUNT      = 0
        self.spinner_char   = ['', '▏', '▎', '▍', '▌', '▋', '▊', '▉', '▊', '▋', '▌', '▍', '▎', '▏',]
        self.COLORS         = ['\033[0;31m', '\033[0;32m', '\033[0;33m']
        self.COLORCOUNT     = 0
        self.COLOR          = self.COLORS[self.COLORCOUNT]
        self.RESET          = '\033[0m'
        self.TSTAGE         = 'file'
        self.TSTART         = False
        self.TSTOP          = False
        self.STARTED        = False
        ## Timer
        # self.START          = time.perf_counter()
        # self.END            = 0
 
    def run(self, main):
        self.output(main)

    def output(self, main):
        def color_change():
            if self.COLORCOUNT == 3:
                self.COLORCOUNT = 0
            self.COLOR = self.COLORS[self.COLORCOUNT]
            self.COLORCOUNT += 1
            return
        
        def spinner(self, main):
            if main.STAGE == 'Finished' and main.PRINTED:
                sys.exit(0)
            if main.ERROR:
                # main.LOG.error_out(main, 'TransRate2', 'TransRate2 failed')
                sys.exit()
            # elif main.STAGE == 'Finished' and not main.PRINTED:
            #     return
            while not main.STAGE:
                pass
            while not main.STARTED:
                pass
            if not self.STARTED:
                self.STAGE = main.STAGE
                self.STARTED = True

            if self.SPINCOUNT == 0:
                color_change()
            self.SPINCOUNT = self.SPINCOUNT + 1 if self.SPINCOUNT < (len(self.spinner_char) - 1) else 0                

            def running():
                out = f'   {self.COLOR}█ {self.RESET}{self.STAGE} Running... {self.COLOR}{self.spinner_char[self.SPINCOUNT]}{self.RESET}'
                print(out, end='\r')

            def finished():
                time_elapsed = round(main.TIMES[self.STAGE], 1)
                stage_string = f'{self.STAGE} Finished'
                out = f'   {main.COLOR}░ {self.RESET}{stage_string:<32}({time_elapsed:.1f}s){" " * 30}'
                print(out)
                self.STAGE = main.STAGE
                main.PRINTED = True

            def transrate():
                if not self.TSTART:
                    print(f'   {main.COLOR}░ {self.RESET}{self.STAGE}{" " * 30}')
                    self.TSTART = True

                elif not main.TSTOP:
                    if self.TSTAGE != main.TSTAGE and main.TSTAGE:
                        time_elapsed = round(main.TIMES[self.TSTAGE], 1)
                        out = f'    {main.COLOR}•  {self.RESET}{self.DESC[self.TSTAGE]:<30}({time_elapsed}s){" " * 30}'
                        print(out)
                        self.TSTAGE = main.TSTAGE
                        main.TPRINTED = True
                    else:
                        desc   = self.DESC[main.TSTAGE]
                        out = f'    {self.COLOR}•  {self.RESET}{desc:<26} {self.COLOR}{self.spinner_char[self.SPINCOUNT]}{self.RESET}'
                        print(out, end='\r')
                        self.TSTAGE = main.TSTAGE
                elif main.TSTAGE:
                    time_elapsed      = round(main.TIMES[self.TSTAGE], 1)
                    stage_string      = f'TransRate2 Finished'
                    print(f'    {main.COLOR}•  {self.RESET}{self.DESC[self.TSTAGE]:<30}({time_elapsed}s){" " * 30}{self.RESET}')
                    print(f'   {main.COLOR}░ {self.RESET}{stage_string:<32}({main.TIMES[main.STAGE]:.1f}s){" " * 30}{self.RESET}')
                    main.TSTAGE = None
                    main.TEND  = True

                    return
                else:
                    self.STAGE = main.STAGE
                    return

            if self.STAGE != 'TransRate2':
                if self.STAGE == main.STAGE:
                    running()
                else:
                    if self.STAGE == 'TransRate2':
                        return
                    else:
                        finished()
            else:
                transrate()

        spinner(self, main)