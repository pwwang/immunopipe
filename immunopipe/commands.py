"""Some other commands for CLI"""
from rich.console import Console
from rich.table import Table

from .utils import text_dedent

def full_opts(params):
    for key, val in params.params.items():
        val.show = True
        if val.type == 'ns':
            for param in val.decendents():
                param.show = True

    params.print_help()
