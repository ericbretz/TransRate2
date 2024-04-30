import os
import sys
import subprocess

def req_check(star = False, bt2 = False, snap=False):
    color          = '\033[0;33m'
    tools          = {
          'salmon'      : 'https://github.com/COMBINE-lab/salmon',
          'samtools'    : 'https://github.com/samtools/samtools',
          'diamond'     : 'https://github.com/bbuchfink/diamond',
        }
    if star:
        tools['STAR'] = 'https://github.com/alexdobin/STAR'
    if bt2:
        tools['bowtie2'] = 'https://github.com/BenLangmead/bowtie2'
    if snap:
        tools['snap-aligner'] = 'https://github.com/amplab/snap'


    missing     = []
    toollabel   = f'{color}  ┌{"─" * 29}\033[m  Requirements  {color}{"─" * 29}┐\033[m'
    bottomlabel = f'{color}  └{"─" * 74}┘\033[m'

    def tool_installed(tool):
        check = subprocess.Popen(['which', tool], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout,stderr = check.communicate()
        return len(stdout.decode('utf-8')) > 0

    for t in tools:
        if not tool_installed(t):
            missing.append(t)

    if len(missing) > 0:
        print(toollabel)
        print(f'{color}  │\033[m  The following tools are missing:{color}{" " * 40}│\033[m')
        for m in missing:
            namelen = len(m)
            linklen = len(tools[m])
            print(f'{color}  │    \033[m  {m}{" " * (15 - namelen)}{tools[m].capitalize()}{" " * (51 - linklen)}{color}  │\033[m')
        print(bottomlabel)
        sys.exit(1)
