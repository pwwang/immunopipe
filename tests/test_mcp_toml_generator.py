"""Tests for MCP TOML generator module."""

from pathlib import Path

import pytest
from immunopipe.mcp.options import OptionsDiscovery
from immunopipe.mcp.toml_generator import (
    _toml_loads,
    ConfigSection,
    TOMLGenerator,
    ConfigTemplateGenerator
)
# Required for test ordering
pytest_order = 100

CONFIGDIR = Path(__file__).parent / "running" / "configs"


class TestConfigSection:
    """Test ConfigSection dataclass."""

    def test_config_section_creation(self):
        """Test creating a ConfigSection instance."""
        section = ConfigSection(
            name="test_section",
            content={"key": "value"},
            description="Test description"
        )

        assert section.name == "test_section"
        assert section.content == {"key": "value"}
        assert section.description == "Test description"


class TestTOMLGenerator:
    """Test TOMLGenerator class."""

    def test_generate_full_config(self):
        """Test generating a complete TOML configuration."""
        generator = TOMLGenerator()

        pipeline_options = {"name": "test_pipeline", "outdir": "./output"}
        process_configs = {
            "SampleInfo": {
                "in": {"infile": ["sample_info.txt"]}
            }
        }
        gbatch_options = {"project": "my-project", "region": "us-central1"}

        config = generator.generate_full_config(
            pipeline_options=pipeline_options,
            process_configs=process_configs,
            gbatch_options=gbatch_options,
            description="Test configuration"
        )

        assert "Test configuration" in config
        assert "name = \"test_pipeline\"" in config
        assert "outdir = \"./output\"" in config
        assert "[SampleInfo]" in config
        assert "[cli-gbatch]" in config
        assert "project = \"my-project\"" in config

    def test_generate_pipeline_section(self):
        """Test generating pipeline options section."""
        generator = TOMLGenerator()

        options = {"name": "test", "forks": 2}
        section = generator.generate_pipeline_section(options)

        assert "Pipeline Options" in section
        assert "name = \"test\"" in section
        assert "forks = 2" in section

    def test_generate_process_section(self):
        """Test generating process configuration section."""
        generator = TOMLGenerator()

        config = {"envs": {"param1": "value1"}}
        section = generator.generate_process_section("TestProcess", config)

        assert "TestProcess process configuration" in section
        assert "[TestProcess]" in section
        assert "[TestProcess.envs]" in section
        assert "param1 = \"value1\"" in section

    def test_generate_gbatch_section(self):
        """Test generating gbatch configuration section."""
        generator = TOMLGenerator()

        options = {"project": "test-project", "region": "us-west1"}
        section = generator.generate_gbatch_section(options)

        assert "Google Batch configuration" in section
        assert "[cli-gbatch]" in section
        assert "project = \"test-project\"" in section
        assert "region = \"us-west1\"" in section

    def test_merge_configs(self):
        """Test merging two configurations."""
        generator = TOMLGenerator()

        base_config = '''
name = "base_pipeline"
outdir = "./base_output"

[SampleInfo]
cache = true
'''

        new_config = '''
name = "new_pipeline"
forks = 4

[SampleInfo]
order = 1

[NewProcess]
enabled = true
'''

        merged = generator.merge_configs(base_config, new_config)
        merged_dict = _toml_loads(merged)

        # Check that new values override base values
        assert merged_dict["name"] == "new_pipeline"
        assert merged_dict["outdir"] == "./base_output"  # preserved
        assert merged_dict["forks"] == 4  # added

        # Check that nested dictionaries are merged
        assert merged_dict["SampleInfo"]["cache"] is True  # preserved
        assert merged_dict["SampleInfo"]["order"] == 1  # added

        # Check that new sections are added
        assert "NewProcess" in merged_dict
        assert merged_dict["NewProcess"]["enabled"] is True

    def test_validate_config(self):
        """Test configuration validation."""
        generator = TOMLGenerator()

        # Valid configuration
        valid_config = '''
name = "test_pipeline"
outdir = "./output"

[SampleInfo]
cache = true

[SampleInfo.envs]
sep = "\\t"
'''

        is_valid, errors = generator.validate_config(valid_config)
        assert is_valid is True
        assert len(errors) == 0

        # Invalid TOML
        invalid_toml = '''
name = "unclosed_string
'''

        is_valid, errors = generator.validate_config(invalid_toml)
        assert is_valid is False
        assert len(errors) > 0
        assert "Validation error" in errors[0]

    def test_deep_merge(self):
        """Test deep merge functionality."""
        generator = TOMLGenerator()

        base = {
            "level1": {
                "level2": {
                    "key1": "value1",
                    "key2": "value2"
                },
                "other_key": "other_value"
            },
            "root_key": "root_value"
        }

        update = {
            "level1": {
                "level2": {
                    "key2": "new_value2",
                    "key3": "value3"
                },
                "new_key": "new_value"
            },
            "new_root_key": "new_root_value"
        }

        result = generator._deep_merge(base, update)

        assert result["level1"]["level2"]["key1"] == "value1"  # preserved
        assert result["level1"]["level2"]["key2"] == "new_value2"  # updated
        assert result["level1"]["level2"]["key3"] == "value3"  # added
        assert result["level1"]["other_key"] == "other_value"  # preserved
        assert result["level1"]["new_key"] == "new_value"  # added
        assert result["root_key"] == "root_value"  # preserved
        assert result["new_root_key"] == "new_root_value"  # added

    def test_extract_section(self):
        """Test extracting a specific section."""
        generator = TOMLGenerator()

        config = '''
name = "test_pipeline"

[SampleInfo]
cache = true

[SampleInfo.envs]
param1 = "value1"

[OtherProcess]
enabled = false
'''

        section = generator.extract_section(config, "SampleInfo")
        section_dict = _toml_loads(section)

        assert "SampleInfo" in section_dict
        assert section_dict["SampleInfo"]["cache"] is True
        assert section_dict["SampleInfo"]["envs"]["param1"] == "value1"

        # Test non-existent section
        missing_section = generator.extract_section(config, "NonExistent")
        assert missing_section is None

    def test_format_config(self):
        """Test configuration formatting."""
        generator = TOMLGenerator()

        config = '''
name="test"
forks=2
[SampleInfo]
cache=true
[cli-gbatch]
project="test-project"
'''

        formatted = generator.format_config(config, add_comments=True)

        assert "Pipeline Configuration Options" in formatted
        assert "SampleInfo Process Configuration" in formatted
        assert "Google Batch Configuration" in formatted

        # Test without comments
        formatted_no_comments = generator.format_config(config, add_comments=False)
        assert "Pipeline Configuration Options" not in formatted_no_comments
        assert "name = \"test\"" in formatted_no_comments


class TestConfigTemplateGenerator:
    """Test ConfigTemplateGenerator class."""

    def test_generate_basic_template(self):
        """Test generating basic template."""
        generator = ConfigTemplateGenerator()
        template = generator.generate_basic_template()

        assert "Basic immunopipe configuration template" in template
        assert "name = \"immunopipe_analysis\"" in template
        assert "outdir = \"./immunopipe_output\"" in template
        assert "[SampleInfo]" in template
        assert "infile" in template and "path/to/sample_info.txt" in template

    def test_generate_tcr_analysis_template(self):
        """Test generating TCR analysis template."""
        generator = ConfigTemplateGenerator()
        template = generator.generate_tcr_analysis_template()

        assert "TCR analysis configuration template" in template
        assert "name = \"tcr_analysis\"" in template
        assert "[TOrBCellSelection]" in template
        assert "[CDR3Clustering]" in template
        assert "[ClonalStats]" in template

    def test_generate_gbatch_template(self):
        """Test generating Google Batch template."""
        generator = ConfigTemplateGenerator()
        template = generator.generate_gbatch_template()

        assert "Google Batch configuration" in template
        assert "[cli-gbatch]" in template
        assert "project = \"your-google-project\"" in template
        assert "region = \"us-central1\"" in template
        assert "machine_type = \"e2-standard-4\"" in template


class TestTCRTemplateOptions:
    """The TCR template must only use the options of the installed pipeline."""

    def test_tcr_template_uses_real_torbcellselection_envs(self):
        """Regression: `cell_type` is not an option of TOrBCellSelection."""
        template = ConfigTemplateGenerator().generate_tcr_analysis_template()
        config = _toml_loads(template)

        envs = config["TOrBCellSelection"]["envs"]
        processes = OptionsDiscovery().get_process_options()
        real_envs = processes["TOrBCellSelection"]["envs"]

        assert envs, "The template should indicate how to select T cells"
        assert "cell_type" not in envs
        assert set(envs) <= set(real_envs), sorted(set(envs) - set(real_envs))

    def test_tcr_template_validates(self):
        """The generated template must pass the server's own validation."""
        template = ConfigTemplateGenerator().generate_tcr_analysis_template()
        is_valid, errors = TOMLGenerator().validate_config(template)

        assert is_valid is True, errors


class TestValidateConfigAgainstPipeline:
    """`validate_config` must check a config against the installed pipeline."""

    def test_invented_sections_and_envs_are_rejected(self):
        config = """
name = "test_pipeline"
outdir = "./output"

[NotAProcess]
envs = {foo = 1}

[TOrBCellSelection.envs]
cell_type = "T"
"""
        is_valid, errors = TOMLGenerator().validate_config(config)

        assert is_valid is False
        assert any("NotAProcess" in error for error in errors), errors
        assert any("cell_type" in error for error in errors), errors

    def test_real_options_at_the_wrong_level_are_warnings(self):
        """pipen silently ignores `min_cells` outside [SeuratPreparing.envs]."""
        config = """
[SeuratPreparing]
min_cells = 3
"""
        is_valid, errors = TOMLGenerator().validate_config(config)

        assert is_valid is True
        assert any(
            error.startswith("Warning:") and "min_cells" in error for error in errors
        ), errors

    def test_real_options_are_accepted(self):
        config = """
name = "test_pipeline"
outdir = "./output"
forks = 2

[SampleInfo.in]
infile = "sample_info.txt"

[TOrBCellSelection.envs]
selector = "Clonotype_Pct > 0.25"
"""
        is_valid, errors = TOMLGenerator().validate_config(config)

        assert is_valid is True, errors

    def test_named_registries_are_accepted(self):
        """`cases`/`stats`-like options take user-chosen names as their keys."""
        config = """
[ClonalStats.envs.cases."A case"]
viz_type = "volume"

[ClonalStats.envs.cases."Another case"]
viz_type = "volume"

[SampleInfo.envs.stats."A stat"]
plot_type = "pie"

[SeuratClusterStats.envs.stats."A cluster stat"]
plot_type = "bar"

[SeuratClusterStats.envs.features."A feature plot"]
plot_type = "violin"
"""
        is_valid, errors = TOMLGenerator().validate_config(config)

        assert is_valid is True, errors

    def test_unknown_pipeline_option_is_rejected(self):
        is_valid, errors = TOMLGenerator().validate_config("bogus_option = 1\n")

        assert is_valid is False
        assert any("bogus_option" in error for error in errors), errors

    def test_repo_configs_have_no_false_positives(self):
        """The repo configs validate, except the genuinely broken ones.

        `malformed.config.toml` is not valid TOML, and the two SeuratPreparing
        configs set `envs.cell_qc_per_sample`, which the installed biopipen
        does not have (pipen silently ignores it and passes it to the job).
        """
        stale = {
            "SeuratPreparing.config.toml",
            "SeuratPreparing_sct.config.toml",
        }

        unexpected = {}
        for path in sorted(CONFIGDIR.glob("*.toml")):
            is_valid, errors = TOMLGenerator().validate_config(path.read_text())

            if path.name == "malformed.config.toml":
                assert is_valid is False, path.name
            elif path.name in stale:
                assert is_valid is False, path.name
                assert all("cell_qc_per_sample" in error for error in errors), (
                    path.name,
                    errors,
                )
            elif not is_valid:
                unexpected[path.name] = errors

        assert unexpected == {}


if __name__ == "__main__":
    pytest.main([__file__])
