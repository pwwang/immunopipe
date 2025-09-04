"""MCP server for immunopipe configuration file composition."""

import sys
import argparse
from typing import List, Optional
import asyncio
import logging

from .server import McpServer

sys.excepthook = sys.__excepthook__


def main(argv: Optional[List[str]] = None) -> None:
    """Main entry point for the MCP server."""
    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser(
        description="MCP server for immunopipe configuration file composition"
    )
    parser.add_argument(
        "--transport",
        choices=["stdio", "http"],
        default="stdio",
        help="Transport method (stdio for VSCode, http for testing)"
    )
    parser.add_argument(
        "--port",
        type=int,
        default=8000,
        help="Port to run the HTTP server on (ignored for stdio)"
    )
    parser.add_argument(
        "--host",
        type=str,
        default="localhost",
        help="Host to run the HTTP server on (ignored for stdio)"
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging level"
    )

    args = parser.parse_args(argv)

    if args.transport == "stdio":
        # For VSCode integration - use stdio transport
        from .mcp_server import main_stdio
        asyncio.run(main_stdio())
    else:
        # For testing - use HTTP transport
        logging.basicConfig(
            level=getattr(logging, args.log_level),
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )

        server = McpServer(host=args.host, port=args.port)
        asyncio.run(server.run())


if __name__ == "__main__":
    main()
