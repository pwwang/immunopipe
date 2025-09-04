# Immunopipe MCP Server

This module provides a Model Context Protocol (MCP) server for composing immunopipe configuration files using natural language descriptions.

## Features

- **Hierarchical Tools**: Tools can call other tools for complex configuration generation
- **Natural Language Processing**: Describe your analysis needs in plain English
- **TOML Generation**: Generate valid TOML configuration files
- **Template Support**: Pre-built templates for common analysis types
- **Configuration Merging**: Combine partial configurations into complete files
- **Process Discovery**: Automatically discover available pipeline processes
- **Validation**: Validate generated configurations

## Available Tools

### Discovery Tools
- `list_pipeline_options` - List all pipeline configuration options
- `list_processes` - List all available immunopipe processes  
- `list_gbatch_options` - List Google Batch configuration options
- `get_process_details` - Get detailed information about a specific process

### Generation Tools
- `generate_pipeline_config` - Generate pipeline-level configuration
- `generate_process_config` - Generate configuration for a specific process
- `generate_gbatch_config` - Generate Google Batch configuration
- `generate_full_config` - Generate a complete configuration file

### Template Tools
- `generate_basic_template` - Generate a basic configuration template
- `generate_tcr_template` - Generate a template for TCR analysis
- `generate_gbatch_template` - Generate a template for Google Batch execution

### Manipulation Tools
- `merge_configs` - Merge two configuration files or strings
- `validate_config` - Validate a TOML configuration
- `format_config` - Format and enhance a TOML configuration

### Assistant Tools
- `help_generate_config` - Get help for generating configuration from natural language
- `suggest_processes` - Suggest processes based on analysis type

## Usage

### Command Line

Start the MCP server:

```bash
python -m immunopipe.mcp --port 8000 --host localhost
```

### Python API

```python
from immunopipe.mcp.server import McpServer, SimpleMCPClient
import asyncio

async def main():
    server = McpServer()
    client = SimpleMCPClient(server)
    
    # Generate a basic template
    result = await client.call_tool("generate_basic_template")
    print(result["result"]["content"])
    
    # Get help for TCR analysis
    result = await client.call_tool("help_generate_config", {
        "description": "I want to analyze single-cell TCR data"
    })
    print(result["result"]["content"]["suggestions"])

asyncio.run(main())
```

### Example Workflow

1. **Describe your analysis**:
   ```python
   result = await client.call_tool("help_generate_config", {
       "description": "I want to perform TCR analysis with clustering"
   })
   ```

2. **Get process suggestions**:
   ```python
   result = await client.call_tool("suggest_processes", {
       "analysis_type": "tcr"
   })
   ```

3. **Generate configuration**:
   ```python
   result = await client.call_tool("generate_full_config", {
       "pipeline_options": {"name": "my_tcr_analysis"},
       "processes": {"TOrBCellSelection": {"envs": {"cell_type": "T"}}},
       "description": "TCR analysis configuration"
   })
   ```

## Configuration Structure

The generated TOML files contain three main sections:

### Pipeline Options
Global settings that control the pipeline behavior:
```toml
name = "my_analysis"
outdir = "./output"
forks = 2
```

### Process Configurations
Settings for individual processes:
```toml
[SampleInfo]
cache = true

[SampleInfo.in]
infile = ["sample_info.txt"]

[SampleInfo.envs]
exclude_cols = "TCRData,BCRData,RNAData"
```

### Google Batch Options
For running on Google Cloud:
```toml
[cli-gbatch]
project = "my-google-project"
region = "us-central1"
machine_type = "e2-standard-4"
```

## Development

### Running Tests

```bash
poetry run pytest tests/test_mcp_*.py -v
```

### Running Examples

```bash
# Automated demonstration
python -m immunopipe.mcp.example

# Interactive mode
python -m immunopipe.mcp.example interactive
```

### Architecture

- `options.py` - Discovers available configuration options from immunopipe
- `toml_generator.py` - Generates and manipulates TOML configuration files
- `tools.py` - Implements the hierarchical MCP tools
- `server.py` - MCP server implementation with JSON-RPC protocol
- `example.py` - Demonstration and testing utilities

## Integration with IDEs

This MCP server can be integrated with IDEs like VSCode through MCP client extensions, enabling natural language configuration generation directly in your development environment.

## Contributing

When adding new tools:

1. Add the tool function to `ImmunopipeConfigTools` class
2. Register it in the `tools` dictionary
3. Add the tool schema to `get_available_tools()`
4. Write tests for the new functionality

## License

This module is part of immunopipe and follows the same license terms.