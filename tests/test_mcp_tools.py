"""Tests for MCP tools module."""
import pytest
from unittest.mock import patch, MagicMock
from immunopipe.mcp.tools import MCPToolResult, ImmunopipeConfigTools
# Required for test ordering
pytest_order = 100


class TestMCPToolResult:
    """Test MCPToolResult dataclass."""

    def test_mcp_tool_result_creation(self):
        """Test creating MCPToolResult instance."""
        result = MCPToolResult(
            success=True,
            content={"key": "value"},
            message="Operation successful",
            tool_calls=[{"tool": "test_tool", "params": {}}],
        )

        assert result.success is True
        assert result.content == {"key": "value"}
        assert result.message == "Operation successful"
        assert len(result.tool_calls) == 1


class TestImmunopipeConfigTools:
    """Test ImmunopipeConfigTools class."""

    def test_initialization(self):
        """Test tools initialization."""
        tools = ImmunopipeConfigTools()

        assert tools.options_discovery is not None
        assert tools.toml_generator is not None
        assert tools.template_generator is not None
        assert len(tools.tools) > 0

        # Check that all expected tools are registered
        expected_tools = [
            "list_pipeline_options",
            "list_processes",
            "list_gbatch_options",
            "get_process_details",
            "generate_pipeline_config",
            "generate_process_config",
            "generate_gbatch_config",
            "generate_full_config",
            "generate_basic_template",
            "generate_tcr_template",
            "generate_gbatch_template",
            "merge_configs",
            "validate_config",
            "format_config",
            "help_generate_config",
            "suggest_processes",
        ]

        for tool in expected_tools:
            assert tool in tools.tools

    def test_get_available_tools(self):
        """Test getting available tools list."""
        tools = ImmunopipeConfigTools()
        available_tools = tools.get_available_tools()

        assert len(available_tools) > 0

        # Check structure of tool definitions
        for tool in available_tools:
            assert "name" in tool
            assert "description" in tool
            assert "inputSchema" in tool
            assert "type" in tool["inputSchema"]
            assert "properties" in tool["inputSchema"]
            assert "required" in tool["inputSchema"]

    @pytest.mark.asyncio
    async def test_execute_tool_unknown(self):
        """Test executing unknown tool."""
        tools = ImmunopipeConfigTools()
        result = await tools.execute_tool("unknown_tool", {})

        assert result.success is False
        assert "Unknown tool" in result.message

    @pytest.mark.asyncio
    @patch.object(ImmunopipeConfigTools, "list_pipeline_options")
    async def test_execute_tool_success(self, mock_tool_func):
        """Test successful tool execution."""
        mock_result = MCPToolResult(success=True, content={"test": "data"})
        mock_tool_func.return_value = mock_result

        tools = ImmunopipeConfigTools()
        result = await tools.execute_tool("list_pipeline_options", {})

        assert result.success is True
        assert result.content == {"test": "data"}
        mock_tool_func.assert_called_once_with()

    @pytest.mark.asyncio
    @patch("immunopipe.mcp.tools.OptionsDiscovery")
    async def test_list_pipeline_options(self, mock_options_discovery):
        """Test list_pipeline_options tool."""
        # Mock the options discovery
        mock_instance = MagicMock()
        mock_option = MagicMock()
        mock_option.description = "Test option"
        mock_option.type = "str"
        mock_option.default = "default"
        mock_option.required = False
        mock_option.choices = None

        mock_instance.get_pipeline_options.return_value = {"test_option": mock_option}
        mock_options_discovery.return_value = mock_instance

        tools = ImmunopipeConfigTools()
        result = await tools.list_pipeline_options()

        assert result.success is True
        assert "test_option" in result.content
        assert result.content["test_option"]["description"] == "Test option"

    @pytest.mark.asyncio
    @patch("immunopipe.mcp.tools.OptionsDiscovery")
    async def test_list_processes(self, mock_options_discovery):
        """Test list_processes tool."""
        mock_instance = MagicMock()
        mock_instance.get_process_options.return_value = {
            "SampleInfo": {
                "description": "Sample information process",
                "envs": {"param1": {"default": "value1"}},
                "requires": [],
            }
        }
        mock_options_discovery.return_value = mock_instance

        tools = ImmunopipeConfigTools()
        result = await tools.list_processes()

        assert result.success is True
        assert "SampleInfo" in result.content
        assert (
            result.content["SampleInfo"]["description"] == "Sample information process"
        )

    @pytest.mark.asyncio
    async def test_generate_basic_template(self):
        """Test generate_basic_template tool."""
        with patch.object(ImmunopipeConfigTools, "__init__", return_value=None):
            tools = ImmunopipeConfigTools()
            tools.template_generator = MagicMock()
            tools.template_generator.generate_basic_template.return_value = (
                "basic template content"
            )

            result = await tools.generate_basic_template()

            assert result.success is True
            assert result.content == "basic template content"
            assert "Generated basic template successfully" in result.message

    @pytest.mark.asyncio
    async def test_merge_configs(self):
        """Test merge_configs tool."""
        with patch.object(ImmunopipeConfigTools, "__init__", return_value=None):
            tools = ImmunopipeConfigTools()
            tools.toml_generator = MagicMock()
            tools.toml_generator.merge_configs.return_value = "merged config content"

            result = await tools.merge_configs("base config", "new config")

            assert result.success is True
            assert result.content == "merged config content"
            assert "Merged configurations successfully" in result.message

    @pytest.mark.asyncio
    async def test_validate_config(self):
        """Test validate_config tool."""
        with patch.object(ImmunopipeConfigTools, "__init__", return_value=None):
            tools = ImmunopipeConfigTools()
            tools.toml_generator = MagicMock()
            tools.toml_generator.validate_config.return_value = (True, [])

            result = await tools.validate_config("config content")

            assert result.success is True
            assert result.content["valid"] is True
            assert result.content["errors"] == []

    @pytest.mark.asyncio
    async def test_help_generate_config_tcr(self):
        """Test help_generate_config with TCR keywords."""
        with patch.object(ImmunopipeConfigTools, "__init__", return_value=None):
            tools = ImmunopipeConfigTools()

            result = await tools.help_generate_config("I want to analyze TCR data")

            assert result.success is True
            assert "suggestions" in result.content
            assert "recommended_tools" in result.content
            assert len(result.content["suggestions"]) > 0
            assert any(
                "TCR" in suggestion for suggestion in result.content["suggestions"]
            )

    @pytest.mark.asyncio
    async def test_help_generate_config_clustering(self):
        """Test help_generate_config with clustering keywords."""
        with patch.object(ImmunopipeConfigTools, "__init__", return_value=None):
            tools = ImmunopipeConfigTools()

            result = await tools.help_generate_config("I need clustering analysis")

            assert result.success is True
            assert "suggestions" in result.content
            assert any(
                "cluster" in suggestion.lower()
                for suggestion in result.content["suggestions"]
            )

    @pytest.mark.asyncio
    async def test_suggest_processes_tcr(self):
        """Test suggest_processes for TCR analysis."""
        with patch.object(ImmunopipeConfigTools, "__init__", return_value=None):
            tools = ImmunopipeConfigTools()

            result = await tools.suggest_processes("tcr")

            assert result.success is True
            assert "suggested_processes" in result.content
            assert "analysis_type" in result.content
            assert result.content["analysis_type"] == "tcr"

            suggested = result.content["suggested_processes"]
            assert "SampleInfo" in suggested
            assert "CDR3Clustering" in suggested
            assert "ClonalStats" in suggested

    @pytest.mark.asyncio
    async def test_suggest_processes_clustering(self):
        """Test suggest_processes for clustering analysis."""
        with patch.object(ImmunopipeConfigTools, "__init__", return_value=None):
            tools = ImmunopipeConfigTools()

            result = await tools.suggest_processes("clustering")

            assert result.success is True
            suggested = result.content["suggested_processes"]
            assert "SeuratClustering" in suggested
            assert "SeuratSubClustering" in suggested
            assert "CellTypeAnnotation" in suggested

    @pytest.mark.asyncio
    async def test_generate_full_config(self):
        """Test generate_full_config tool."""
        with patch.object(ImmunopipeConfigTools, "__init__", return_value=None):
            tools = ImmunopipeConfigTools()
            tools.toml_generator = MagicMock()
            tools.toml_generator.generate_full_config.return_value = (
                "full config content"
            )

            pipeline_opts = {"name": "test"}
            processes = {"SampleInfo": {}}
            gbatch_opts = {"project": "test-project"}

            result = await tools.generate_full_config(
                pipeline_options=pipeline_opts,
                processes=processes,
                gbatch_options=gbatch_opts,
                description="Test config",
            )

            assert result.success is True
            assert result.content == "full config content"
            tools.toml_generator.generate_full_config.assert_called_once_with(
                pipeline_options=pipeline_opts,
                process_configs=processes,
                gbatch_options=gbatch_opts,
                description="Test config",
            )

    @pytest.mark.asyncio
    async def test_get_process_details_existing(self):
        """Test get_process_details for existing process."""
        with patch.object(ImmunopipeConfigTools, "__init__", return_value=None):
            tools = ImmunopipeConfigTools()
            tools.options_discovery = MagicMock()
            tools.options_discovery.get_process_options.return_value = {
                "SampleInfo": {
                    "description": "Sample info process",
                    "envs": {},
                    "requires": [],
                }
            }

            result = await tools.get_process_details("SampleInfo")

            assert result.success is True
            assert result.content["description"] == "Sample info process"

    @pytest.mark.asyncio
    async def test_get_process_details_not_found(self):
        """Test get_process_details for non-existent process."""
        with patch.object(ImmunopipeConfigTools, "__init__", return_value=None):
            tools = ImmunopipeConfigTools()
            tools.options_discovery = MagicMock()
            tools.options_discovery.get_process_options.return_value = {}

            result = await tools.get_process_details("NonExistentProcess")

            assert result.success is False
            assert "not found" in result.message


if __name__ == "__main__":
    pytest.main([__file__])
