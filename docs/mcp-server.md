# MCP Server

The **Immunopipe MCP Server** is a Model Context Protocol (MCP) server that provides intelligent configuration generation for immunopipe pipelines. It allows AI assistants like Claude to understand and generate complex TOML configuration files for single-cell RNA-seq and TCR/BCR analysis workflows.

## Overview

The MCP server exposes a comprehensive set of tools for:

- **Discovery**: List available pipeline options, processes, and Google Batch configurations
- **Generation**: Create TOML configuration files, templates, and process-specific configs
- **Manipulation**: Merge, validate, and format configuration files
- **AI-Assisted Configuration**: Generate configurations from natural language descriptions

## Installation

The MCP server is included with immunopipe. You can run it directly using:

```bash
# Run with stdio transport (for VSCode/Claude integration)
immunopipe mcp --transport stdio

# Run with HTTP transport (for testing)
immunopipe mcp --transport http --port 8000
```

## VSCode Integration

To integrate with VSCode and Claude, add the following to your VSCode settings:

```json
{
  "mcp.servers": {
    "immunopipe": {
      "command": "python",
      "args": ["-m", "immunopipe.mcp", "--transport", "stdio"],
      "env": {}
    }
  }
}
```

## Available Tools

### Discovery Tools

#### `list_pipeline_options`
Lists all available pipeline-level configuration options.

**Parameters**: None

**Returns**: Dictionary of pipeline options with descriptions, types, defaults, and requirements.

#### `list_processes`
Lists all available immunopipe processes.

**Parameters**: None

**Returns**: Dictionary of processes with descriptions, environment variables, and requirements.

#### `list_gbatch_options`
Lists all Google Batch configuration options.

**Parameters**: None

**Returns**: Dictionary of Google Batch options with descriptions, types, and defaults.

#### `get_process_details`
Gets detailed information about a specific process.

**Parameters**:
- `process_name` (string): Name of the process to get details for

**Returns**: Detailed process information including configuration options.

### Generation Tools

#### `generate_full_config`
Generates a complete immunopipe TOML configuration file.

**Parameters**:
- `pipeline_options` (object, optional): Pipeline-level options
- `processes` (object, optional): Process configurations
- `gbatch_options` (object, optional): Google Batch options
- `description` (string, optional): Description for the configuration file

**Returns**: Complete TOML configuration as a string.

#### `generate_pipeline_config`
Generates pipeline-level configuration section.

**Parameters**:
- `options` (object): Pipeline options to include

**Returns**: Pipeline configuration TOML section.

#### `generate_process_config`
Generates configuration for a specific process.

**Parameters**:
- `process_name` (string): Name of the process
- `config` (object): Process configuration options

**Returns**: Process configuration TOML section.

#### `generate_gbatch_config`
Generates Google Batch configuration section.

**Parameters**:
- `options` (object): Google Batch options

**Returns**: Google Batch configuration TOML section.

### Template Tools

#### `generate_basic_template`
Generates a basic immunopipe configuration template.

**Parameters**: None

**Returns**: Basic TOML template suitable for standard scRNA-seq analysis.

#### `generate_tcr_template`
Generates a template optimized for TCR analysis.

**Parameters**: None

**Returns**: TOML template configured for single-cell TCR/BCR analysis.

#### `generate_gbatch_template`
Generates a template for Google Batch execution.

**Parameters**: None

**Returns**: TOML template configured for cloud execution via Google Batch.

### Manipulation Tools

#### `merge_configs`
Merges two TOML configuration files or strings.

**Parameters**:
- `base_config` (string): Base configuration (TOML string)
- `new_config` (string): New configuration to merge (TOML string)

**Returns**: Merged TOML configuration.

#### `validate_config`
Validates a TOML configuration file.

**Parameters**:
- `config_content` (string): TOML configuration content to validate

**Returns**: Validation results with success status and error messages.

#### `format_config`
Formats and enhances a TOML configuration.

**Parameters**:
- `config_content` (string): TOML configuration content to format
- `add_comments` (boolean, optional): Whether to add helpful comments

**Returns**: Formatted TOML configuration with optional comments.

### AI-Assisted Tools

#### `configure_immunopipe`
Generates immunopipe configuration from natural language requirements. This is the most intelligent tool that can understand complex analysis requirements and automatically select appropriate processes.

**Parameters**:
- `requirements` (string): Natural language description of analysis needs
- `specific_processes` (array, optional): Specific process names to configure
- `parameters` (object, optional): Specific parameter values to set

**Examples**:
- "clustering with resolution 0.5"
- "TCR analysis with clustering"
- "set SeuratPreparing QC parameters"

**Returns**: Complete TOML configuration tailored to the requirements.

#### `help_generate_config`
Provides configuration suggestions based on natural language descriptions.

**Parameters**:
- `description` (string): Natural language description of analysis goals

**Returns**: Suggestions and recommended tools for the analysis.

#### `suggest_processes`
Suggests processes based on analysis type.

**Parameters**:
- `analysis_type` (string): Type of analysis (e.g., 'tcr', 'bcr', 'clustering', 'differential')

**Returns**: List of recommended processes for the analysis type.

## Usage Examples

### Basic Usage

```python
from immunopipe.mcp.tools import ImmunopipeConfigTools

# Initialize the tools
tools = ImmunopipeConfigTools()

# Generate a basic template
result = await tools.generate_basic_template()
print(result.content)
```

### Natural Language Configuration

```python
# Configure immunopipe using natural language
result = await tools.configure_immunopipe(
    requirements="I want to analyze single-cell TCR data with clustering resolution 0.8"
)
print(result.content)  # Complete TOML configuration
```

### Process Discovery

```python
# List available processes
processes = await tools.list_processes()
print(f"Found {len(processes.content)} processes")

# Get details for a specific process
details = await tools.get_process_details("SeuratClustering")
print(details.content)
```

### Configuration Manipulation

```python
# Generate base config
base = await tools.generate_basic_template()

# Generate TCR-specific additions
tcr_config = await tools.generate_tcr_template()

# Merge configurations
merged = await tools.merge_configs(base.content, tcr_config.content)
print(merged.content)
```

## Integration with Claude/AI Assistants

The MCP server is designed to work seamlessly with AI assistants. With `immunopipe` installed,
you can add the MCP server to your Claude or VSCode setup to enable intelligent configuration generation:

```shell
claude mcp add --transport http immunopipe-mcp http://localhost:80000
```

```json
{
  "mcp": {
    "servers": {
        "immunopipe": {
            "command": "python",
            "args": ["-m", "immunopipe", "mcp", "--transport", "stdio"]
        }
    }
}
```

Here are some example prompts that work well:

### Configuration Generation
- "Generate an immunopipe configuration for TCR analysis"
- "Create a config for single-cell clustering with resolution 0.5"
- "Set up immunopipe for differential expression analysis"

### Process Selection
- "What processes do I need for BCR analysis?"
- "Suggest processes for clustering analysis"
- "Show me options for Google Batch execution"

### Configuration Management
- "Validate this immunopipe configuration"
- "Format and add comments to my config file"
- "Merge these two immunopipe configurations"

## Architecture

The MCP server consists of several key components:

- **`McpServer`**: Main server implementation with JSON-RPC over HTTP
- **`MCPServer`**: Standards-compliant MCP server for VSCode integration
- **`ImmunopipeConfigTools`**: Core tool collection with hierarchical capabilities
- **`OptionsDiscovery`**: Discovers available configuration options and processes
- **`TOMLGenerator`**: Generates and manipulates TOML configuration files
- **`ProcessDocumentationExtractor`**: Extracts process documentation for intelligent configuration

## Error Handling

The MCP server provides comprehensive error handling:

- **Tool Validation**: Validates tool parameters and provides helpful error messages
- **Configuration Validation**: Checks TOML syntax and immunopipe-specific requirements
- **Process Discovery**: Handles missing or unavailable processes gracefully
- **Natural Language Processing**: Provides fallback configurations when analysis fails

## Development and Testing

To test the MCP server locally:

```bash
# Run the example demonstration
python -m immunopipe.mcp.example

# Run interactive demo
python -m immunopipe.mcp.example interactive

# Test with HTTP transport
immunopipe mcp --transport http --port 8000
```

The server includes comprehensive test coverage and example usage patterns to help with development and integration.

## Performance

The MCP server is optimized for:

- **Fast Discovery**: Efficient caching of process and option information
- **Incremental Configuration**: Tools can be chained for complex configurations
- **Memory Efficiency**: Lazy loading of process documentation
- **Scalability**: Async/await support for concurrent tool execution

## Support and Troubleshooting

Common issues and solutions:

1. **Process Not Found**: Ensure immunopipe is properly installed and processes are available
2. **Configuration Validation Errors**: Check TOML syntax and parameter validity
3. **Natural Language Analysis Failures**: Use more specific descriptions or explicit process names
4. **VSCode Integration Issues**: Verify MCP server configuration in VSCode settings

For additional help, open an issue on the [GitHub repository](https://github.com/pwwang/immunopipe).
