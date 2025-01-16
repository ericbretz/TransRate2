from rich.console import Console
from rich.theme import Theme


class TerminalRemake:
    def __init__(self):
        self.console = Console()
        return
    
    def main(self):
        self.console.print('Hello, World!', style='bold red')
        return
    
if __name__ == '__main__':
    terminal = TerminalRemake()
    terminal.main()
