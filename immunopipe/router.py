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
        print("Commands:")
        print("  gbatch    Run entire pipeline as a job on Google Cloud Batch")
        print(
            "  mcp       "
            "Run MCP (model context protocol) server for configuration generation"
        )
        print("  utils     Run utility functions ")
        print(
            "            (e.g. checking gene symbols, Seurat object dimensions, etc.)"
        )
        print("")
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
    elif "-v" in sys.argv or "--version" in sys.argv:
        print(__version__)
    else:
        from .pipeline import main as pipeline_main

        pipeline_main()
