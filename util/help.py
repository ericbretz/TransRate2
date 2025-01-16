class HelpText:
    def __init__(self, version):
        self.version = version

    def printLogo(self):
        logo = f"""
 ██████┐██████┐  █████┐ ███┐  ██┐ ██████┐██████┐  █████┐ ██████┐███████┐██████┐
 └─██┌─┘██┌──██┐██┌──██┐████┐ ██│██┌────┘██┌──██┐██┌──██┐└─██┌─┘██┌────┘└────██┐
   ▓▓│  ▓▓▓▓▓▓┌┘▓▓▓▓▓▓▓│▓▓┌▓▓┐▓▓│└▓▓▓▓▓┐ ▓▓▓▓▓▓┌┘▓▓▓▓▓▓▓│  ▓▓│  ▓▓▓▓▓┐    ▓▓▓┌─┘
   ▒▒│  ▒▒┌──▒▒┐▒▒┌──▒▒│▒▒│└▒▒▒▒│ └───▒▒┐▒▒┌──▒▒┐▒▒┌──▒▒│  ▒▒│  ▒▒┌──┘  ▒▒┌──┘
   ░░│  ░░│  ░░│░░│  ░░│░░│ └░░░│░░░░░░┌┘░░│  ░░│░░│  ░░│  ░░│  ░░░░░░░┐░░░░░░░┐
   └─┘  └─┘  └─┘└─┘  └─┘└─┘  └──┘└─────┘ └─┘  └─┘└─┘  └─┘  └─┘  └──────┘└──────┘
             Quality analysis for de-novo transcriptome assemblies
             ░▓▓▓▓▓▓▓◠▓▓▓▓▓▓▓░ ░▓▓▓▓▓▓▓◠▓▓▓▓▓▓▓░ ░▓▓▓▓▓▓▓◠▓▓▓▓▓▓▓░      EC Bretz
                                                                          v{self.version}"""
        print(logo)

    def printBox(self):
        body = """
  ┌─────────────────────────────  Help Options  ─────────────────────────────┐
  │                                                                          │
  ├┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄    Assembly    ┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┤
  │                                                                          │
  │ --assembly         -a        Path to assembly file (FASTA)               │
  │ --left             -l        Path to left reads file (FASTQ)             │
  │ --right            -r        Path to right reads file (FASTQ)            │
  │ --reference        -f        Path to reference file (FASTA)              │
  │ --outdir           -o        Path to output directory                    │
  │ --threads          -t        Number of threads to use                    │
  │                                                                          │
  ├┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄    Aligners    ┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┤
  │                                                                          │
  │ --hisat2           -s        Use HISAT2 aligner                          │
  │ --bowtie2          -b        Use Bowtie2 aligner                         │
  │                                                                          │
  ├┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄     Others     ┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┤
  │                                                                          │
  │ --clutter          -c        Remove intermediate files                   │
  │ --quiet            -q        Supress terminal output                     │
  │ --help             -h        Display this help message                   │
  │ --version          -v        Display version                             │
  │                                                                          │
  └──────────────────────────────────────────────────────────────────────────┘
  ┌─────────────────────────────   Mode Types   ─────────────────────────────┐
  │                                                                          │
  │ -a                  Run assembly analysis only.                          │
  │ -a -l -r            Run assembly with paired-end reads analysis          │
  │ -a -l -r -f         Run assembly with paired-end reads and reference     │
  │ -a -l               Run assembly with single-end reads analysis          │
  │ -a -f               Run assembly with reference analysis                 │
  │                                                                          │
  └──────────────────────────────────────────────────────────────────────────┘
  """
        print(body)

    def run(self):
        self.printLogo()
        self.printBox()
