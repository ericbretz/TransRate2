class Logo:
    def __init__(self, version):
        self.version = version
        self.author = 'EC Bretz'
        self.colors = {
            'red'   : '\033[0;31m',
            'green' : '\033[0;32m',
            'yellow': '\033[0;33m',
            'blue'  : '\033[0;34m',
            'purple': '\033[0;35m',
            'cyan'  : '\033[0;36m',
            'white' : '\033[0;37m',
            'reset' : '\033[0m'
        }
        self.term_colors = [self.colors['yellow'], self.colors['red'], self.colors['green']]
        self.term_reset = self.colors['reset']
        self.logo = self.transrate_logo()

    def transrate_logo(self):
        W = self.colors['reset']
        G = self.colors['green']
        Y = self.colors['yellow']
        R = self.colors['red']

        logo = f'''
 {G}██████{W}┐{G}██████{W}┐  {G}█████{W}┐ {G}███{W}┐  {G}██{W}┐ {G}██████{W}┐{Y}██████{W}┐  {Y}█████{W}┐ {Y}██████{W}┐{Y}███████{W}┐{R}██████{W}┐ 
 └─{G}██{W}┌─┘{G}██{W}┌──{G}██{W}┐{G}██{W}┌──{G}██{W}┐{G}████{W}┐ {G}██{W}│{G}██{W}┌────┘{Y}██{W}┌──{Y}██{W}┐{Y}██{W}┌──{Y}██{W}┐└─{Y}██{W}┌─┘{Y}██{W}┌────┘└────{R}██{W}┐
   {G}▓▓{W}│  {G}▓▓▓▓▓▓{W}┌┘{G}▓▓▓▓▓▓▓{W}│{G}▓▓{W}┌{G}▓▓{W}┐{G}▓▓{W}│└{G}▓▓▓▓▓{W}┐ {Y}▓▓▓▓▓▓{W}┌┘{Y}▓▓▓▓▓▓▓{W}│  {Y}▓▓{W}│  {Y}▓▓▓▓▓{W}┐    {R}▓▓▓{W}┌─┘
   {G}▒▒{W}│  {G}▒▒{W}┌──{G}▒▒{W}┐{G}▒▒{W}┌──{G}▒▒{W}│{G}▒▒{W}│└{G}▒▒▒▒{W}│ └───{G}▒▒{W}┐{Y}▒▒{W}┌──{Y}▒▒{W}┐{Y}▒▒{W}┌──{Y}▒▒{W}│  {Y}▒▒{W}│  {Y}▒▒{W}┌──┘  {R}▒▒{W}┌──┘  
   {G}░░{W}│  {G}░░{W}│  {G}░░{W}│{G}░░{W}│  {G}░░{W}│{G}░░{W}│ └{G}░░░{W}│{G}░░░░░░{W}┌┘{Y}░░{W}│  {Y}░░{W}│{Y}░░{W}│  {Y}░░{W}│  {Y}░░{W}│  {Y}░░░░░░░{W}┐{R}░░░░░░░{W}┐
 {W}  └─┘  └─┘  └─┘└─┘  └─┘└─┘  └──┘└─────┘ └─┘  └─┘└─┘  └─┘  └─┘  └──────┘└──────┘
{W}{"Quality analysis for de-novo transcriptome assemblies":^80}
{W}             {G}░▓▓▓▓▓▓▓◠▓▓▓▓▓▓▓░ {Y}░▓▓▓▓▓▓▓◠▓▓▓▓▓▓▓░ {R}░▓▓▓▓▓▓▓◠▓▓▓▓▓▓▓░ {W}   EC Bretz
{W}{self.version:>78}\033[0m
        '''
        return logo
    