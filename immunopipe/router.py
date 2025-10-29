"""Routing the program execution based on the 2nd argument

(with the first being the program name).

- If the 2nd argument is "gbatch", we are employing pipen-cli-gbatch to run the
whole pipeline on Google Cloud Batch.
- If the 2nd argument is "mcp", we are running a MCP (model context protocol) server.
- Otherwise, we are running the pipeline directly with given configuration file or
arguments.
"""

import sys
import asyncio
from .version import desc, __version__

COMMANDS = {
    "gbatch": "Run entire pipeline as a job on Google Cloud Batch",
    "mcp": "Run MCP (model context protocol) server for configuration generation",
    "utils": (
        "Run utility functions",
        "(e.g. checking gene symbols, Seurat object dimensions, etc.)",
    ),
    "help": "Show the help message for specific processes",
}


def _help_commands() -> str:
    """Generate help message for commands.

    Returns:
        The help message string.
    """
    help_msg = "Commands:\n"
    for cmd, help_line in COMMANDS.items():
        if isinstance(help_line, tuple):
            for i, line in enumerate(help_line):
                if i == 0:
                    help_msg += f"  {cmd:10} {line}\n"
                else:
                    help_msg += f"             {line}\n"
        else:
            help_msg += f"  {cmd:10} {help_line}\n"
    return help_msg


def run():
    """Route the program execution based on the 2nd argument."""
    if (
        len(sys.argv) < 2
        or (
            len(sys.argv) == 2
            and "-h" not in sys.argv
            and "--help" not in sys.argv
            and "-h+" not in sys.argv
            and "--help+" not in sys.argv
            and "-v" not in sys.argv
            and "--version" not in sys.argv
            and sys.argv[1] not in COMMANDS
            and not sys.argv[1].startswith("@")
        )
    ):
        print(
            "\033[1;4mUsage\033[0m: immunopipe "
            "[command|options to run pipeline directly]"
        )
        print(
            "       immunopipe <-h|--help>  # Show options to run pipeline directly"
        )
        print("")
        print(desc)
        print("")
        print(_help_commands())
        sys.exit(0)

    command = sys.argv[1]

    if command == "gbatch":
        from .gbatch import main as gbatch_main
        asyncio.run(gbatch_main(sys.argv[2:]))

    elif command == "mcp":
        from .mcp import main as mcp_main
        mcp_main(sys.argv[2:])

    elif command == "utils":
        from .cli_utils import main as utils_main
        utils_main(sys.argv[2:])

    elif command == "help":
        from .cli_help import main as help_main
        help_main(sys.argv[2:])

    elif "-v" in sys.argv or "--version" in sys.argv:
        print(__version__)

    else:
        from .pipeline import main as pipeline_main
        pipeline_main()
