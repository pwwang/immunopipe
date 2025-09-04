"""MCP Server implementation for immunopipe configuration."""

import json
import asyncio
import logging
from typing import Dict, Any, List

from .tools import ImmunopipeConfigTools

logger = logging.getLogger(__name__)


class MCPRequest:
    """MCP request wrapper."""

    def __init__(
        self, method: str, params: Dict[str, Any] = None, request_id: str = None
    ):
        self.method = method
        self.params = params or {}
        self.request_id = request_id


class MCPResponse:
    """MCP response wrapper."""

    def __init__(
        self, result: Any = None, error: Dict[str, Any] = None, request_id: str = None
    ):
        self.result = result
        self.error = error
        self.request_id = request_id

    def to_dict(self) -> Dict[str, Any]:
        """Convert response to dictionary."""
        response = {}

        if self.request_id is not None:
            response["id"] = self.request_id

        if self.error:
            response["error"] = self.error
        else:
            response["result"] = self.result

        return response


class McpServer:
    """MCP Server for immunopipe configuration generation."""

    def __init__(self, host: str = "localhost", port: int = 8000):
        self.host = host
        self.port = port
        self.tools = ImmunopipeConfigTools()
        self.server_info = {
            "name": "immunopipe-config-mcp",
            "version": "1.0.0",
            "description": "MCP server for immunopipe configuration file composition",
            "author": "Immunopipe Team",
            "capabilities": {
                "tools": True,
                "hierarchical_tools": True,
                "configuration_generation": True,
                "toml_manipulation": True,
            },
        }

    async def run(self):
        """Run the MCP server."""
        logger.info(f"Starting MCP server on {self.host}:{self.port}")

        # For now, we'll implement a simple JSON-RPC over HTTP
        # In a real implementation, you'd use a proper MCP transport
        try:
            await self._start_server()
        except KeyboardInterrupt:
            logger.info("Server stopped by user")
        except Exception as e:
            logger.error(f"Server error: {e}")

    async def _start_server(self):
        """Start the server (placeholder for actual MCP transport)."""
        # This is a simplified implementation
        # In practice, you would use the official MCP transport libraries

        logger.info("MCP Server running...")
        logger.info(f"Available tools: {list(self.tools.tools.keys())}")

        # Keep server running
        while True:
            await asyncio.sleep(1)

    async def handle_request(self, request_data: str) -> str:
        """Handle incoming MCP request."""
        try:
            request_dict = json.loads(request_data)
            request = MCPRequest(
                method=request_dict.get("method"),
                params=request_dict.get("params", {}),
                request_id=request_dict.get("id"),
            )

            response = await self._process_request(request)
            return json.dumps(response.to_dict())

        except Exception as e:
            error_response = MCPResponse(
                error={"code": -32603, "message": f"Internal error: {str(e)}"}
            )
            return json.dumps(error_response.to_dict())

    async def _process_request(self, request: MCPRequest) -> MCPResponse:
        """Process MCP request and return response."""

        if request.method == "initialize":
            return MCPResponse(result=self.server_info, request_id=request.request_id)

        elif request.method == "tools/list":
            tools = self.tools.get_available_tools()
            return MCPResponse(result={"tools": tools}, request_id=request.request_id)

        elif request.method == "tools/call":
            return await self._handle_tool_call(request)

        else:
            return MCPResponse(
                error={
                    "code": -32601,
                    "message": f"Method not found: {request.method}",
                },
                request_id=request.request_id,
            )

    async def _handle_tool_call(self, request: MCPRequest) -> MCPResponse:
        """Handle tool call request."""
        try:
            tool_name = request.params.get("name")
            arguments = request.params.get("arguments", {})

            if not tool_name:
                return MCPResponse(
                    error={"code": -32602, "message": "Missing tool name"},
                    request_id=request.request_id,
                )

            # Execute the tool
            result = await self.tools.execute_tool(tool_name, arguments)

            if result.success:
                response_content = {
                    "content": result.content,
                    "message": result.message,
                }

                # Handle hierarchical tool calls
                if result.tool_calls:
                    response_content["tool_calls"] = result.tool_calls

                return MCPResponse(
                    result=response_content, request_id=request.request_id
                )
            else:
                return MCPResponse(
                    error={"code": -32603, "message": result.message},
                    request_id=request.request_id,
                )

        except Exception as e:
            return MCPResponse(
                error={"code": -32603, "message": f"Tool execution error: {str(e)}"},
                request_id=request.request_id,
            )

    def get_server_info(self) -> Dict[str, Any]:
        """Get server information."""
        return self.server_info


class SimpleMCPClient:
    """Simple MCP client for testing the server."""

    def __init__(self, server: McpServer):
        self.server = server

    async def call_tool(
        self, tool_name: str, arguments: Dict[str, Any] = None
    ) -> Dict[str, Any]:
        """Call a tool on the server."""
        request_data = {
            "method": "tools/call",
            "params": {"name": tool_name, "arguments": arguments or {}},
            "id": "test-request",
        }

        response_str = await self.server.handle_request(json.dumps(request_data))
        response = json.loads(response_str)

        return response

    async def list_tools(self) -> List[Dict[str, Any]]:
        """List available tools."""
        request_data = {"method": "tools/list", "params": {}, "id": "list-tools"}

        response_str = await self.server.handle_request(json.dumps(request_data))
        response = json.loads(response_str)

        if "result" in response:
            return response["result"].get("tools", [])
        else:
            return []

    async def initialize(self) -> Dict[str, Any]:
        """Initialize connection with server."""
        request_data = {"method": "initialize", "params": {}, "id": "init"}

        response_str = await self.server.handle_request(json.dumps(request_data))
        response = json.loads(response_str)

        return response


# Example usage and testing functions
async def test_server():
    """Test the MCP server functionality."""
    server = McpServer()
    client = SimpleMCPClient(server)

    # Initialize
    init_response = await client.initialize()
    print("Initialize response:", init_response)

    # List tools
    tools = await client.list_tools()
    print(f"Available tools: {len(tools)}")
    for tool in tools[:3]:  # Show first 3 tools
        print(f"  - {tool['name']}: {tool['description']}")

    # Test a simple tool
    result = await client.call_tool("generate_basic_template")
    print("Basic template result:", result.get("result", {}).get("message", ""))

    # Test pipeline options
    result = await client.call_tool("list_pipeline_options")
    if result.get("result", {}).get("content"):
        options_count = len(result["result"]["content"])
        print(f"Found {options_count} pipeline options")

    return server


if __name__ == "__main__":
    # Test the server
    async def main():
        await test_server()

    asyncio.run(main())
