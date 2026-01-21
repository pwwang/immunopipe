# IMMUNOPIPE/MCP/ PROJECT KNOWLEDGE BASE

**Generated:** 2026-01-20
**Module:** MCP Server for AI Configuration Generation

## OVERVIEW

MCP server enabling natural language → TOML config generation via hierarchical tools (2937 lines, 8 files).

## WHERE TO LOOK

| Task | Location | Notes |
|------|----------|-------|
| Tool organization | `tools.py` (1107 lines) | ImmunopipeConfigTools class, tool registry dict |
| Add new MCP tool | `tools.py` line 34 | Add to `self.tools` dict + implement async method |
| Tool schema | `tools.py` line 59 | `get_available_tools()` returns JSON Schema |
| Stdio transport (VSCode) | `mcp_server.py` line 40 | `run_stdio()` for IDE integration |
| HTTP transport (testing) | `__init__.py` line 60 | Server mode via argparse |
| Config discovery | `options.py` (486 lines) | Pipen args + process auto-discovery |
| Process env extraction | `doc_extractor.py` (295 lines) | pipen-annotate integration |
| TOML generation | `toml_generator.py` (322 lines) | Diot/simpleconf utilities |

## CONVENTIONS

### Tool Architecture Pattern
```python
# tools.py - Centralized tool registry
self.tools = {
    "tool_name": self.tool_method,  # Direct mapping
}

async def execute_tool(self, tool_name, params):
    tool_func = self.tools[tool_name]
    result = await tool_func(**params)
    return MCPToolResult(success=True, content=result)
```

### Adding New MCP Tools
1. Implement async method in `ImmunopipeConfigTools` class
2. Add entry to `self.tools` dict (line 34)
3. Define JSON Schema in `get_available_tools()` (line 59)
4. Return `MCPToolResult(success=bool, content=Any, message=str)`

### Hierarchical Tool Calls
- High-level tools (e.g., `configure_immunopipe`) call lower-level tools internally
- Result includes `tool_calls` list for chaining: `MCPToolResult(..., tool_calls=[...])`
- Example: `configure_immunopipe` → `list_processes` → `generate_process_config`

### Transport Modes
| Mode | Entry Point | Use Case |
|------|-------------|----------|
| stdio | `mcp_server.py:main_stdio()` | VSCode/Claude Desktop (stdio transport) |
| http | `server.py:McpServer.run()` | Testing/dev (HTTP on port 8000) |

### Config Discovery Mechanism
- **Pipeline options**: `pipen_args.defaults.PIPEN_ARGS` introspection
- **Processes**: Temporarily reload `immunopipe.processes` with mock config to expose conditional processes
- **Process envs**: `pipen-annotate` extracts structured metadata from process classes

### TOML Generation Pattern
```python
# toml_generator.py - Diot for TOML serialization
from diot import Diot
toml_str = Diot(config_dict).to_toml()

# Section merging uses simpleconf.Config.load()
```

## ANTI-PATTERNS

**None unique to this module** - Follows MCP protocol + project conventions.
