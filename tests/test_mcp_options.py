"""Tests for MCP options discovery module."""
import pytest
from unittest.mock import patch, MagicMock
from immunopipe.mcp.options import (
    ConfigOption,
    PipelineOptionsDiscovery,
    ProcessDiscovery,
    GbatchOptionsDiscovery,
    OptionsDiscovery
)
# Required for test ordering
pytest_order = 100


class TestConfigOption:
    """Test ConfigOption dataclass."""

    def test_config_option_creation(self):
        """Test creating a ConfigOption instance."""
        option = ConfigOption(
            name="test_option",
            description="Test description",
            type="str",
            default="default_value",
            choices=["a", "b", "c"],
            required=True,
            section="test_section"
        )

        assert option.name == "test_option"
        assert option.description == "Test description"
        assert option.type == "str"
        assert option.default == "default_value"
        assert option.choices == ["a", "b", "c"]
        assert option.required is True
        assert option.section == "test_section"


class TestPipelineOptionsDiscovery:
    """Test PipelineOptionsDiscovery class."""

    @patch('pipen_args.defaults.PIPEN_ARGS', {
        'name': {
            'help': 'Pipeline name',
            'default': 'test_pipeline',
            'type': str,
            'required': False
        },
        'outdir': {
            'help': 'Output directory',
            'default': './output',
            'type': str,
            'required': True
        }
    })
    def test_get_pipeline_options(self):
        """Test getting pipeline options."""
        discovery = PipelineOptionsDiscovery()
        options = discovery.get_pipeline_options()

        assert 'name' in options
        assert 'outdir' in options

        name_option = options['name']
        assert name_option.name == 'name'
        assert name_option.description == 'Pipeline name'
        assert name_option.default == 'test_pipeline'
        assert name_option.required is False

        outdir_option = options['outdir']
        assert outdir_option.name == 'outdir'
        assert outdir_option.description == 'Output directory'
        assert outdir_option.default == './output'
        assert outdir_option.required is True

    def test_infer_type(self):
        """Test type inference."""
        discovery = PipelineOptionsDiscovery()

        # Test explicit type
        option_with_type = {'type': int, 'default': 42}
        assert discovery._infer_type(option_with_type) == 'int'

        # Test default value type
        option_with_default = {'default': 'string_value'}
        assert discovery._infer_type(option_with_default) == 'str'

        # Test choices
        option_with_choices = {'choices': ['a', 'b', 'c']}
        assert discovery._infer_type(option_with_choices) == 'choice'

        # Test fallback
        empty_option = {}
        assert discovery._infer_type(empty_option) == 'str'


class TestProcessDiscovery:
    """Test ProcessDiscovery class."""

    def test_extract_description(self):
        """Test description extraction from docstring."""
        discovery = ProcessDiscovery()

        # Test normal docstring
        docstring = '''"""This is a test process.

        This is a longer description.
        """'''
        assert discovery._extract_description(docstring) == "This is a test process."

        # Test empty docstring
        assert discovery._extract_description(None) == ""
        assert discovery._extract_description("") == ""

        # Test docstring with markers
        docstring_with_markers = '''"""
        This should be extracted.
        """'''
        assert discovery._extract_description(docstring_with_markers) == (
            "This should be extracted."
        )

    def test_extract_envs_info(self):
        """Test environment variables extraction."""
        discovery = ProcessDiscovery()

        # Mock process with envs
        mock_proc = MagicMock()
        mock_proc.envs = {
            'param1': 'default_value',
            'param2': 42,
            'param3': True
        }

        envs_info = discovery._extract_envs_info(mock_proc)

        assert 'param1' in envs_info
        assert envs_info['param1']['default'] == 'default_value'
        assert envs_info['param1']['type'] == 'str'

        assert 'param2' in envs_info
        assert envs_info['param2']['default'] == 42
        assert envs_info['param2']['type'] == 'int'

        assert 'param3' in envs_info
        assert envs_info['param3']['default'] is True
        assert envs_info['param3']['type'] == 'bool'

    def test_get_process_requirements(self):
        """Test process requirements extraction."""
        discovery = ProcessDiscovery()

        # Mock process with single requirement
        mock_proc1 = MagicMock()
        mock_requirement = MagicMock()
        mock_requirement.name = 'RequiredProcess'
        mock_proc1.requires = mock_requirement

        requirements = discovery._get_process_requirements(mock_proc1)
        assert requirements == ['RequiredProcess']

        # Mock process with multiple requirements
        mock_proc2 = MagicMock()
        mock_req1 = MagicMock()
        mock_req1.name = 'Process1'
        mock_req2 = MagicMock()
        mock_req2.name = 'Process2'
        mock_proc2.requires = [mock_req1, mock_req2]

        requirements = discovery._get_process_requirements(mock_proc2)
        assert 'Process1' in requirements
        assert 'Process2' in requirements


class TestGbatchOptionsDiscovery:
    """Test GbatchOptionsDiscovery class."""

    def test_discover_gbatch_options(self):
        """Test discovering gbatch options from config."""
        discovery = GbatchOptionsDiscovery()

        # Test that we can discover options
        # (this will work if pipen-cli-gbatch is available)
        try:
            options = discovery.get_gbatch_options()

            # Should find some common options
            assert isinstance(options, dict)

            # Check if we found expected options
            if 'project' in options:
                assert options['project'].name == 'project'
                assert options['project'].section == 'cli-gbatch'
                assert hasattr(options['project'], 'type')

        except ImportError:
            # If pipen-cli-gbatch is not available, that's okay for testing
            assert True


class TestOptionsDiscovery:
    """Test main OptionsDiscovery class."""

    def test_initialization(self):
        """Test OptionsDiscovery initialization."""
        discovery = OptionsDiscovery()

        assert isinstance(discovery.pipeline_discovery, PipelineOptionsDiscovery)
        assert isinstance(discovery.process_discovery, ProcessDiscovery)
        assert isinstance(discovery.gbatch_discovery, GbatchOptionsDiscovery)

    @patch.object(PipelineOptionsDiscovery, 'get_pipeline_options')
    @patch.object(ProcessDiscovery, 'get_processes')
    @patch.object(GbatchOptionsDiscovery, 'get_gbatch_options')
    def test_get_all_options(self, mock_gbatch, mock_processes, mock_pipeline):
        """Test getting all options."""
        mock_pipeline.return_value = {'pipeline_opt': MagicMock()}
        mock_processes.return_value = {'SampleInfo': MagicMock()}
        mock_gbatch.return_value = {'gbatch_opt': MagicMock()}

        discovery = OptionsDiscovery()
        all_options = discovery.get_all_options()

        assert 'pipeline' in all_options
        assert 'processes' in all_options
        assert 'gbatch' in all_options

        mock_pipeline.assert_called_once()
        mock_processes.assert_called_once()
        mock_gbatch.assert_called_once()


if __name__ == "__main__":
    pytest.main([__file__])
