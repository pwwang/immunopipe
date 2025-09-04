"""Proper MCP server implementation for VSCode integration."""

import json
import sys
import asyncio
import logging
from typing import Dict, Any, Optional

from .tools import ImmunopipeConfigTools

logger = logging.getLogger(__name__)


class MCPServer:
    """Standards-compliant MCP server for VSCode integration."""

    def __init__(self):
        self.tools = ImmunopipeConfigTools()
        self.server_info = {
            "name": "immunopipe-mcp",
            "version": "1.0.0",
            "protocolVersion": "2024-11-05",
            "serverInfo": {
                "name": "immunopipe-mcp",
                "version": "1.0.0",
                "description": (
                    "Immunopipe configuration generator for scRNA-seq and "
                    "scTCR/BCR analysis"
                ),
            },
            "capabilities": {
                "tools": {"listChanged": False},
                "logging": {"level": "info"},
                "completion": {"supportsCompletion": True},
                "experimental": {},
                "resources": {},
            },
        }

    async def run_stdio(self):
        """Run MCP server over stdio for VSCode integration."""
        # Server is ready, no need to log startup message

        try:
            while True:
                # Read JSON-RPC message from stdin
                line = await asyncio.get_event_loop().run_in_executor(
                    None, sys.stdin.readline
                )

                if not line:
                    break

                line = line.strip()
                if not line:
                    continue

                try:
                    request = json.loads(line)
                    response = await self._handle_request(request)

                    if response:
                        # Write response to stdout
                        print(json.dumps(response), flush=True)

                except json.JSONDecodeError as e:
                    logger.error(f"Invalid JSON received: {e}")
                    error_response = {
                        "jsonrpc": "2.0",
                        "id": None,
                        "error": {"code": -32700, "message": "Parse error"},
                    }
                    print(json.dumps(error_response), flush=True)
                except Exception as e:
                    logger.error(f"Error handling request: {e}")

        except Exception as e:
            logger.error(f"Server error: {e}")

    async def _handle_request(
        self, request: Dict[str, Any]
    ) -> Optional[Dict[str, Any]]:
        """Handle MCP request and return response."""
        method = request.get("method")
        params = request.get("params", {})
        request_id = request.get("id")

        logger.debug(f"Handling request: {method}")

        try:
            if method == "initialize":
                return {"jsonrpc": "2.0", "id": request_id, "result": self.server_info}

            elif method == "initialized":
                # No response needed for initialized notification
                return None

            elif method == "tools/list":
                tools = self.tools.get_available_tools()
                return {"jsonrpc": "2.0", "id": request_id, "result": {"tools": tools}}

            elif method == "tools/call":
                tool_name = params.get("name")
                arguments = params.get("arguments", {})

                if not tool_name:
                    return {
                        "jsonrpc": "2.0",
                        "id": request_id,
                        "error": {"code": -32602, "message": "Missing tool name"},
                    }

                # Execute the tool
                result = await self.tools.execute_tool(tool_name, arguments)

                if result.success:
                    response_content = [
                        {
                            "type": "text",
                            "text": (
                                result.content
                                if isinstance(result.content, str)
                                else json.dumps(result.content, indent=2)
                            ),
                        }
                    ]

                    return {
                        "jsonrpc": "2.0",
                        "id": request_id,
                        "result": {"content": response_content, "isError": False},
                    }
                else:
                    return {
                        "jsonrpc": "2.0",
                        "id": request_id,
                        "error": {"code": -32603, "message": result.message},
                    }

            elif method == "ping":
                return {"jsonrpc": "2.0", "id": request_id, "result": {}}

            elif method == "logging/setLevel":
                # Handle VSCode's logging level setting
                level = params.get("level", "info")
                logger.info(f"Setting log level to: {level}")
                return {"jsonrpc": "2.0", "id": request_id, "result": {}}

            elif method == "notifications/cancelled":
                # Handle cancellation notifications
                return None

            elif method == "completion/complete":
                # Handle completion requests
                return {
                    "jsonrpc": "2.0",
                    "id": request_id,
                    "result": {
                        "completion": {"values": [], "total": 0, "hasMore": False}
                    },
                }

            else:
                return {
                    "jsonrpc": "2.0",
                    "id": request_id,
                    "error": {"code": -32601, "message": f"Method not found: {method}"},
                }

        except Exception as e:
            logger.error(f"Error processing request {method}: {e}")
            return {
                "jsonrpc": "2.0",
                "id": request_id,
                "error": {"code": -32603, "message": f"Internal error: {str(e)}"},
            }


async def main_stdio():
    """Main entry point for stdio MCP server."""
    # Set up minimal logging to stderr so it doesn't interfere with stdout
    logging.basicConfig(
        level=logging.WARNING,  # Reduce logging verbosity
        format="%(name)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stderr)],
    )

    # Only log important startup info
    logger = logging.getLogger(__name__)
    logger.warning("Immunopipe MCP server starting...")

    server = MCPServer()
    await server.run_stdio()


if __name__ == "__main__":
    asyncio.run(main_stdio())
