"""TOML configuration file generation and manipulation."""

from typing import Dict, Any, List, Optional, Tuple
from pathlib import Path
from dataclasses import dataclass
import logging

from diot import Diot
from simpleconf import Config

from .options import get_config_surface

logger = logging.getLogger(__name__)

# Prefix of the messages that report a problem that does not invalidate a
# configuration (e.g. an option used at the wrong level)
WARNING_PREFIX = "Warning: "


def _toml_dumps(data: Dict[str, Any]) -> str:
    """Dump dictionary to TOML string with proper formatting."""
    return Diot(data).to_toml()


def _toml_loads(toml_str: str) -> Dict[str, Any]:
    """Load TOML string into a dictionary."""
    return Config.load(toml_str, loader="tomls")


@dataclass
class ConfigSection:
    """Represents a configuration section."""

    name: str
    content: Dict[str, Any]
    description: Optional[str] = None


class TOMLGenerator:
    """Generate and manipulate TOML configuration files for immunopipe."""

    # Options of a process section, besides `envs`/`in`/`out`. These are the
    # options of `pipen.Proc`, which are the ones pipen reads from a process
    # section of a configuration file.
    PROCESS_OPTIONS = frozenset(
        (
            "name",
            "desc",
            "envs_depth",
            "cache",
            "dirsig",
            "export",
            "error_strategy",
            "num_retries",
            "template",
            "template_opts",
            "forks",
            "input",
            "input_data",
            "lang",
            "order",
            "output",
            "output_flatten",
            "plugin_opts",
            "requires",
            "scheduler",
            "scheduler_opts",
            "script",
            "submission_batch",
            "workdir",
            # Not an option of pipen, but used by the config generator of
            # this server to materialize an otherwise empty section
            "enabled",
        )
    )

    def __init__(self):
        pass

    def generate_full_config(
        self,
        pipeline_options: Dict[str, Any] = None,
        process_configs: Dict[str, Dict[str, Any]] = None,
        gbatch_options: Dict[str, Any] = None,
        description: str = "Immunopipe configuration file",
    ) -> str:
        """Generate a complete TOML configuration file."""

        config = {}
        lines = [f"# {description}\n"]

        # Add pipeline options
        if pipeline_options:
            lines.append("# Pipeline Options")
            for key, value in pipeline_options.items():
                if value is not None:
                    config[key] = value

        # Add process configurations
        if process_configs:
            for process_name, process_config in process_configs.items():
                lines.append(f"\n# {process_name} process configuration")
                # Include process even if config is empty (this enables it in pipeline)
                process_dict = process_config if process_config else {}

                # For processes with only nested config (like 'in', 'out', 'envs'),
                # add an empty process section to ensure [ProcessName] appears
                if process_dict and all(
                    key in ["in", "out", "envs"] for key in process_dict.keys()
                ):
                    # Add an empty value to force creation of [ProcessName] section
                    process_dict["enabled"] = True

                config[process_name] = process_dict

        # Add gbatch options
        if gbatch_options:
            lines.append("\n# Google Batch configuration")
            config["cli-gbatch"] = gbatch_options

        # Generate TOML string
        toml_content = _toml_dumps(config)

        # Combine custom comments with TOML content
        return "\n".join(lines) + "\n\n" + toml_content

    def generate_pipeline_section(self, options: Dict[str, Any]) -> str:
        """Generate only the pipeline options section."""
        if not options:
            return "# No pipeline options specified\n"

        return "# Pipeline Options\n" + _toml_dumps(options)

    def generate_process_section(
        self, process_name: str, config: Dict[str, Any]
    ) -> str:
        """Generate a process configuration section."""
        if not config:
            return f"# [{process_name}]\n# No configuration specified\n"

        result = f"# {process_name} process configuration\n"

        # Ensure we create a [ProcessName] section header when needed
        config_copy = dict(config)
        if config_copy and all(
            key in ["in", "out", "envs"] for key in config_copy.keys()
        ):
            # Add an empty value to force creation of [ProcessName] section
            config_copy["enabled"] = True

        process_config = {process_name: config_copy}
        result += _toml_dumps(process_config)

        return result

    def generate_gbatch_section(self, options: Dict[str, Any]) -> str:
        """Generate the cli-gbatch section."""
        if not options:
            return "# [cli-gbatch]\n# No gbatch options specified\n"

        return "# Google Batch configuration\n" + _toml_dumps({"cli-gbatch": options})

    def merge_configs(self, base_config: str, new_config: str) -> str:
        """Merge a new configuration into an existing base configuration."""
        try:
            # Parse both configs
            if base_config.strip():
                base_dict = _toml_loads(base_config)
            else:
                base_dict = {}

            new_dict = _toml_loads(new_config)

            # Merge configurations
            merged = self._deep_merge(base_dict, new_dict)

            return _toml_dumps(merged)

        except Exception as e:
            logger.error(f"Failed to merge configs: {e}")
            return base_config

    def merge_config_files(self, base_file: str, new_config: str) -> str:
        """Merge new configuration into an existing config file."""
        try:
            if Path(base_file).exists():
                with open(base_file, "r") as f:
                    base_content = f.read()
            else:
                base_content = ""

            return self.merge_configs(base_content, new_config)

        except Exception as e:
            logger.error(f"Failed to merge config file {base_file}: {e}")
            return new_config

    def validate_config(self, config_content: str) -> Tuple[bool, List[str]]:
        """Validate TOML configuration content against the installed pipeline.

        A section must be a process of the installed pipeline (or a
        pipeline-level option), and an `envs` key must be an environment
        variable of that process. Options that do exist, but are used at the
        wrong level (e.g. an environment variable set directly in the process
        section, which pipen silently ignores) are reported as warnings,
        prefixed with `Warning:`; they don't make a config invalid.

        Returns:
            A tuple with whether the config is valid, and the messages.
        """
        try:
            config_dict = _toml_loads(config_content)
        except Exception as e:
            return False, [f"Validation error: {e}"]

        # Basic validation
        if not isinstance(config_dict, dict):
            return False, ["Configuration must be a dictionary"]

        surface = get_config_surface()
        errors = []

        for key, value in config_dict.items():
            if key == "cli-gbatch":
                continue  # Options for running on Google Batch

            if key in surface["pipeline_options"]:
                continue  # Pipeline options take values of any type

            if not isinstance(value, dict):
                errors.append(f"Unknown pipeline option: {key}")
                continue

            if key not in surface["processes"]:
                errors.append(f"Unknown process section: [{key}]")
                continue

            errors.extend(
                self._validate_process_section(key, value, surface["processes"][key])
            )

        # Warnings are reported along with the errors, but don't make the
        # config invalid
        is_valid = not [
            error for error in errors if not error.startswith(WARNING_PREFIX)
        ]
        return is_valid, errors

    def _validate_process_section(
        self, process_name: str, section: Dict[str, Any], surface: Dict[str, Any]
    ) -> List[str]:
        """Validate a process section against the options of the process."""
        errors = []
        envs = surface.get("envs") or set()

        for key, value in section.items():
            if key == "envs":
                if not isinstance(value, dict):
                    errors.append(f"Process {process_name}: envs must be a dictionary")
                    continue

                # Nothing below the first level is checked: the keys of a
                # registry-like option (e.g. `ClonalStats.envs.cases`) are
                # user-chosen names
                for env_key in value:
                    if env_key not in envs:
                        errors.append(
                            "Unknown environment variable of process "
                            f"{process_name}: {env_key}"
                        )
            elif key in ("in", "out"):
                # The names of the input/output variables of the process
                known = surface.get("inputs" if key == "in" else "outputs")
                if not known or not isinstance(value, dict):
                    continue

                for io_key in value:
                    if io_key not in known:
                        errors.append(
                            f"Unknown {key} variable of process "
                            f"{process_name}: {io_key}"
                        )
            elif key in self.PROCESS_OPTIONS:
                continue
            elif key in envs:
                errors.append(
                    f"{WARNING_PREFIX}[{process_name}] {key} is an environment "
                    f"variable of process {process_name} and will be ignored, "
                    f"set it as [{process_name}.envs] {key} instead"
                )
            else:
                errors.append(f"Unknown option of process {process_name}: {key}")

        return errors

    def _deep_merge(
        self, base: Dict[str, Any], update: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Deep merge two dictionaries."""
        result = base.copy()

        for key, value in update.items():
            if (
                key in result
                and isinstance(result[key], dict)
                and isinstance(value, dict)
            ):
                result[key] = self._deep_merge(result[key], value)
            else:
                result[key] = value

        return result

    def extract_section(self, config_content: str, section_name: str) -> Optional[str]:
        """Extract a specific section from a TOML config."""
        try:
            config_dict = _toml_loads(config_content)

            if section_name in config_dict:
                section_config = {section_name: config_dict[section_name]}
                return _toml_dumps(section_config)

            return None

        except Exception as e:
            logger.error(f"Failed to extract section {section_name}: {e}")
            return None

    def format_config(self, config_content: str, add_comments: bool = True) -> str:
        """Format and enhance TOML configuration with comments."""
        try:
            config_dict = _toml_loads(config_content)

            if not add_comments:
                return _toml_dumps(config_dict)

            # Add helpful comments
            lines = []

            # Pipeline options section
            pipeline_keys = set(config_dict.keys()) - {"cli-gbatch"}
            process_keys = {
                k for k in pipeline_keys if isinstance(config_dict.get(k), dict)
            }
            simple_keys = pipeline_keys - process_keys

            if simple_keys:
                lines.append("# Pipeline Configuration Options")
                for key in sorted(simple_keys):
                    # Format individual values properly
                    value = config_dict[key]
                    if isinstance(value, str):
                        formatted_value = f'"{value}"'
                    elif isinstance(value, bool):
                        formatted_value = str(value).lower()
                    else:
                        formatted_value = str(value)
                    lines.append(f"{key} = {formatted_value}")
                lines.append("")

            # Process sections
            for process_name in sorted(process_keys):
                lines.append(f"# {process_name} Process Configuration")
                process_config = {process_name: config_dict[process_name]}
                lines.append(_toml_dumps(process_config).strip())
                lines.append("")

            # Gbatch section
            if "cli-gbatch" in config_dict:
                lines.append("# Google Batch Configuration")
                gbatch_config = {"cli-gbatch": config_dict["cli-gbatch"]}
                lines.append(_toml_dumps(gbatch_config).strip())

            return "\n".join(lines)

        except Exception as e:
            logger.error(f"Failed to format config: {e}")
            return config_content


class ConfigTemplateGenerator:
    """Generate configuration templates for common use cases."""

    def __init__(self, toml_generator: TOMLGenerator = None):
        self.toml_gen = toml_generator or TOMLGenerator()

    def generate_basic_template(self) -> str:
        """Generate a basic configuration template."""
        config = {
            "name": "immunopipe_analysis",
            "outdir": "./immunopipe_output",
            "forks": 1,
            "SampleInfo": {"in": {"infile": ["path/to/sample_info.txt"]}},
        }

        return self.toml_gen.generate_full_config(
            pipeline_options={
                "name": config["name"],
                "outdir": config["outdir"],
                "forks": config["forks"],
            },
            process_configs={"SampleInfo": config["SampleInfo"]},
            description="Basic immunopipe configuration template",
        )

    def generate_tcr_analysis_template(self) -> str:
        """Generate a template for TCR analysis."""
        config = {
            "name": "tcr_analysis",
            "outdir": "./tcr_output",
            "TOrBCellSelection": {
                # `selector` is an expression to indicate which cells are
                # T cells, see `immunopipe help TOrBCellSelection`
                "envs": {"selector": "Clonotype_Pct > 0.25"}
            },
            "CDR3Clustering": {},
            "ClonalStats": {},
        }

        return self.toml_gen.generate_full_config(
            pipeline_options={"name": config["name"], "outdir": config["outdir"]},
            process_configs={
                k: v for k, v in config.items() if k not in ["name", "outdir"]
            },
            description="TCR analysis configuration template",
        )

    def generate_gbatch_template(self) -> str:
        """Generate a template for Google Batch execution."""
        gbatch_config = {
            "project": "your-google-project",
            "region": "us-central1",
            "machine_type": "e2-standard-4",
            "disk_size": "100GB",
        }

        return self.toml_gen.generate_gbatch_section(gbatch_config)
