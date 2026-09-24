"""Options discovery for immunopipe configuration."""

import sys
from typing import Dict, Any, List, Optional, Set
from dataclasses import dataclass
import logging

from pipen import Proc
from pipen.defaults import CONFIG as PIPEN_CONFIG
from pipen.utils import LOADING_ARGV0

from .doc_extractor import ProcessDocumentationExtractor

logger = logging.getLogger(__name__)


@dataclass
class ConfigOption:
    """Represents a configuration option."""

    name: str
    description: str
    type: str
    default: Any
    choices: Optional[List[str]] = None
    required: bool = False
    section: Optional[str] = None


class PipelineOptionsDiscovery:
    """Discover pipeline options from pipen_args.defaults."""

    def __init__(self):
        self._pipeline_options: Optional[Dict[str, ConfigOption]] = None

    def get_pipeline_options(self) -> Dict[str, ConfigOption]:
        """Get all pipeline options from pipen_args.defaults."""
        if self._pipeline_options is None:
            self._discover_pipeline_options()
        return self._pipeline_options

    def _discover_pipeline_options(self) -> None:
        """Discover pipeline options from pipen_args.defaults.PIPEN_ARGS."""
        try:
            from pipen_args.defaults import PIPEN_ARGS

            self._pipeline_options = {}

            for name, option in PIPEN_ARGS.items():
                # Extract option details
                config_option = ConfigOption(
                    name=name,
                    description=option.get("help", ""),
                    type=self._infer_type(option),
                    default=option.get("default"),
                    choices=option.get("choices"),
                    required=option.get("required", False),
                    section="pipeline",
                )
                self._pipeline_options[name] = config_option

        except Exception as e:
            logger.error(f"Failed to discover pipeline options: {e}")
            self._pipeline_options = {}

    def _infer_type(self, option: Dict[str, Any]) -> str:
        """Infer the type of an option.

        The types in `PIPEN_ARGS` are not necessarily type objects, they can
        also be strings (e.g. `"auto"`) or callables (e.g. `str.upper`).
        """
        if "type" in option:
            option_type = option["type"]
            # `str(...)` as a fallback for types that have no `__name__`
            return getattr(option_type, "__name__", None) or str(option_type)

        default = option.get("default")
        if default is not None:
            return type(default).__name__

        if "choices" in option:
            return "choice"

        return "str"


class ProcessDiscovery:
    """Discover available processes and their options."""

    def __init__(self):
        self._processes: Optional[Dict[str, Dict[str, Any]]] = None
        self._config_surface: Optional[Dict[str, Dict[str, Set[str]]]] = None

    def get_processes(self) -> Dict[str, Dict[str, Any]]:
        """Get all available processes and their information."""
        if self._processes is None:
            self._discover_processes()
        return self._processes

    def _discover_processes(self) -> None:
        """Discover available processes from the immunopipe pipeline."""
        try:
            # Temporarily set LOADING_ARGV0 to notify immunopipe we're only
            # loading structure
            original_argv0 = sys.argv[0]
            sys.argv[0] = LOADING_ARGV0

            self._processes = {}

            # First, reload the processes module to force all conditional processes to
            # be evaluated
            # We'll do this by temporarily modifying the config to include conditional
            # triggers
            try:
                # Create a mock configuration that enables conditional processes
                import tempfile
                import os

                mock_config_content = """
[TOrBCellSelection]
envs = {}

[ClusterMarkersOfAllCells]
envs = {}

[TopExpressingGenesOfAllCells]
envs = {}

[LoadingRNAFromSeurat.envs]
prepared = false
clustered = false
"""

                # Write mock config to temp file
                with tempfile.NamedTemporaryFile(
                    mode="w", suffix=".toml", delete=False
                ) as f:
                    f.write(mock_config_content)
                    temp_config_file = f.name

                # Temporarily set the config file in sys.argv
                original_argv = sys.argv[:]
                sys.argv = ["immunopipe", f"@{temp_config_file}"]

                try:
                    # Force reimport of processes module with mock config
                    import importlib
                    import immunopipe.processes

                    importlib.reload(immunopipe.processes)

                    # Now get the pipeline with conditional processes loaded
                    from immunopipe.pipeline import Immunopipe

                    pipe = Immunopipe()
                    pipe.build_proc_relationships()

                    # Get all processes from pipeline
                    for proc in pipe.procs:
                        proc_info = {
                            "name": proc.name,
                            "description": self._extract_description(proc.__doc__),
                            "envs": self._extract_envs_info(proc),
                            "requires": self._get_process_requirements(proc),
                            "doc": proc.__doc__ or "",
                        }
                        self._processes[proc.name] = proc_info

                finally:
                    # Restore original argv and clean up
                    sys.argv = original_argv
                    try:
                        os.unlink(temp_config_file)
                    except Exception:
                        pass

            except Exception as e:
                logger.warning(
                    f"Could not load conditional processes with mock config: {e}"
                )

                # Fallback to normal pipeline loading
                from immunopipe.pipeline import Immunopipe

                pipe = Immunopipe()
                pipe.build_proc_relationships()

                for proc in pipe.procs:
                    proc_info = {
                        "name": proc.name,
                        "description": self._extract_description(proc.__doc__),
                        "envs": self._extract_envs_info(proc),
                        "requires": self._get_process_requirements(proc),
                        "doc": proc.__doc__ or "",
                    }
                    self._processes[proc.name] = proc_info

            # Also scan the processes module for any additional process classes
            try:
                import immunopipe.processes

                for attr_name in dir(immunopipe.processes):
                    if (
                        attr_name.startswith("_")
                        or attr_name in self._processes
                        or not hasattr(
                            getattr(immunopipe.processes, attr_name), "__doc__"
                        )
                    ):
                        continue

                    attr = getattr(immunopipe.processes, attr_name)

                    # Check if this is a process class (has process-like attributes)
                    if (
                        hasattr(attr, "envs")
                        or hasattr(attr, "requires")
                        or ("Seurat" in attr_name and "Clustering" in attr_name)
                        or ("Cluster" in attr_name and "Markers" in attr_name)
                        or ("TopExpressing" in attr_name and "Genes" in attr_name)
                    ):

                        proc_info = {
                            "name": attr_name,
                            "description": self._extract_description(attr.__doc__),
                            "envs": (
                                self._extract_envs_info(attr)
                                if hasattr(attr, "envs")
                                else {}
                            ),
                            "requires": (
                                self._get_process_requirements(attr)
                                if hasattr(attr, "requires")
                                else []
                            ),
                            "doc": attr.__doc__ or "",
                        }
                        self._processes[attr_name] = proc_info
                        logger.debug(
                            f"Added process class from module scan: {attr_name}"
                        )

            except Exception as e:
                logger.warning(f"Error scanning processes module: {e}")

            # Restore original argv[0]
            sys.argv[0] = original_argv0

            logger.info(f"Total processes discovered: {len(self._processes)}")

            # Check for specific conditional processes
            conditional_found = []
            conditional_missing = []
            for proc_name in [
                "SeuratClusteringOfAllCells",
                "ClusterMarkersOfAllCells",
                "TopExpressingGenesOfAllCells",
            ]:
                if proc_name in self._processes:
                    conditional_found.append(proc_name)
                else:
                    conditional_missing.append(proc_name)

            if conditional_found:
                logger.info(
                    f"✓ Conditional processes found: {', '.join(conditional_found)}"
                )
            if conditional_missing:
                logger.warning(
                    f"✗ Conditional processes missing: {', '.join(conditional_missing)}"
                )

        except Exception as e:
            logger.error(f"Failed to discover processes: {e}")
            self._processes = {}

    def get_config_surface(self) -> Dict[str, Dict[str, Set[str]]]:
        """Get the real configuration surface of the installed pipeline.

        Returns:
            A mapping from process name to the names of the options that the
            process accepts, so that configurations can be checked against
            what the installed pipeline actually reads:

                {
                    "SampleInfo": {
                        "envs": {"stats", "sep", ...},
                        "inputs": {"infile"},
                        "outputs": {"outfile", "outdir"},
                    },
                    ...
                }
        """
        if self._config_surface is None:
            self._discover_config_surface()
        return self._config_surface

    def _discover_config_surface(self) -> None:
        """Discover the configuration surface from the pipeline and its docs.

        Unlike `_discover_processes()`, the surface is read from a pipeline
        with all of the conditional processes defined, so that any process
        that can be enabled from a configuration file is included.
        """
        self._config_surface = {}

        # Temporarily set LOADING_ARGV0 to notify immunopipe that we are
        # only loading the structure, so that all of the conditional
        # processes are defined
        original_argv = sys.argv[:]
        sys.argv = [LOADING_ARGV0]

        try:
            import importlib
            import immunopipe.processes
            import immunopipe.pipeline

            # The processes are defined conditionally at import time, so
            # both modules need to be reloaded for all of them to show up.
            # `pipeline` is reloaded as well, because it snapshots
            # `processes.start_processes`.
            importlib.reload(immunopipe.processes)
            importlib.reload(immunopipe.pipeline)

            extractor = ProcessDocumentationExtractor()

            try:
                pipe = immunopipe.pipeline.Immunopipe()
                pipe.build_proc_relationships()
                for proc in pipe.procs:
                    self._add_proc_to_surface(proc.__class__, proc.name, extractor)
            except Exception as e:
                logger.warning(f"Could not load all processes from the pipeline: {e}")

            # Also scan the processes module for any additional process
            # classes (e.g. the ones not triggered by the current pipeline).
            # Only the classes defined by the module are considered, to skip
            # the base classes imported from `biopipen`.
            for attr in vars(immunopipe.processes).values():
                if (
                    isinstance(attr, type)
                    and getattr(attr, "__module__", None)
                    == immunopipe.processes.__name__
                    and issubclass(attr, Proc)
                ):
                    self._add_proc_to_surface(attr, None, extractor)

        except Exception as e:
            logger.error(f"Failed to discover configuration surface: {e}")
        finally:
            sys.argv = original_argv

    def _add_proc_to_surface(
        self,
        proc: Any,
        name: Optional[str],
        extractor: ProcessDocumentationExtractor,
    ) -> None:
        """Add the options of a process to the configuration surface."""
        name = name or proc.__name__
        surface = self._config_surface.setdefault(
            name, {"envs": set(), "inputs": set(), "outputs": set()}
        )

        # The defaults defined by the process class ...
        surface["envs"].update(getattr(proc, "envs", None) or {})

        # ... and the ones documented by the process
        annotations = extractor.get_process_annotations(proc)
        surface["envs"].update(annotations.get("Envs", None) or {})
        surface["inputs"].update(annotations.get("Input", None) or {})
        surface["outputs"].update(annotations.get("Output", None) or {})

    def _extract_description(self, docstring: Optional[str]) -> str:
        """Extract a brief description from a docstring."""
        if not docstring:
            return ""

        # Find the first line that's not empty or a docstring marker
        lines = docstring.strip().split("\n")
        for line in lines:
            cleaned = line.strip()
            if not cleaned:
                continue

            # Handle docstring markers at start or end of line
            if cleaned.startswith('"""') or cleaned.startswith("'''"):
                # Check if there's content after the marker on the same line
                for marker in ['"""', "'''"]:
                    if cleaned.startswith(marker):
                        content_after_marker = cleaned[len(marker) :].strip()
                        if content_after_marker and not content_after_marker.endswith(
                            marker
                        ):
                            return content_after_marker
            elif not cleaned.endswith('"""') and not cleaned.endswith("'''"):
                return cleaned

        return ""

    def _extract_envs_info(self, proc) -> Dict[str, Any]:
        """Extract environment variables information from a process."""
        envs_info = {}

        if hasattr(proc, "envs") and proc.envs:
            for key, value in proc.envs.items():
                envs_info[key] = {
                    "default": value,
                    "type": type(value).__name__,
                    "description": f"Environment variable for {key}",
                }

        return envs_info

    def _get_process_requirements(self, proc) -> List[str]:
        """Get the requirements for a process."""
        requirements = []

        if hasattr(proc, "requires"):
            if isinstance(proc.requires, (list, tuple)):
                requirements = [
                    req.name for req in proc.requires if hasattr(req, "name")
                ]
            elif hasattr(proc.requires, "name"):
                requirements = [proc.requires.name]

        return requirements


class GbatchOptionsDiscovery:
    """Discover cli-gbatch options from the gbatch parser."""

    def __init__(self):
        self._gbatch_options: Optional[Dict[str, ConfigOption]] = None

    def get_gbatch_options(self) -> Dict[str, ConfigOption]:
        """Get cli-gbatch options from the gbatch parser."""
        if self._gbatch_options is None:
            self._discover_gbatch_options()
        return self._gbatch_options

    def _discover_gbatch_options(self) -> None:
        """Discover gbatch options from the parser configuration."""
        try:
            # Import gbatch components
            from pathlib import Path
            from simpleconf import Config
            from pipen_cli_gbatch import __file__ as cli_gbatch_file

            # Load cli-gbatch arguments configuration
            cli_gbatch_arg_file = Path(cli_gbatch_file).parent / "daemon_args.toml"
            cli_gbatch_args = Config.load_one(cli_gbatch_arg_file, loader="toml")

            self._gbatch_options = {}

            # Process action options
            if (
                hasattr(cli_gbatch_args, "mutually_exclusive_groups")
                and cli_gbatch_args.mutually_exclusive_groups
            ):
                action_options = cli_gbatch_args.mutually_exclusive_groups[0].arguments
                for opt in action_options:
                    option = self._create_config_option_from_arg(opt, section="Actions")
                    if option:
                        self._gbatch_options[option.name] = option

            # Process regular arguments
            cli_gbatch_arguments = []
            if hasattr(cli_gbatch_args, "arguments"):
                cli_gbatch_arguments.extend(cli_gbatch_args.arguments)
            if hasattr(cli_gbatch_args, "groups"):
                for group in cli_gbatch_args.groups:
                    if hasattr(group, "arguments"):
                        cli_gbatch_arguments.extend(
                            [
                                arg
                                for arg in group.arguments
                                if arg.get("flags") != ["command"]
                            ]
                        )

            for arg in cli_gbatch_arguments:
                # Skip certain options that are handled specially
                arg_flags = arg.get("flags", [])
                if any(
                    flag
                    in [
                        "--name",
                        "--plain",
                        "--workdir",
                        "--jobname-prefix",
                        "--cwd",
                        "--entrypoint",
                        "--commands",
                    ]
                    for flag in arg_flags
                ):
                    continue

                option = self._create_config_option_from_arg(arg, section="cli-gbatch")
                if option:
                    self._gbatch_options[option.name] = option

        except Exception as e:
            logger.error(f"Failed to discover gbatch options: {e}")
            self._gbatch_options = {}

    def _create_config_option_from_arg(
        self, arg_dict: dict, section: str
    ) -> Optional[ConfigOption]:
        """Create a ConfigOption from an argument dictionary."""
        try:
            flags = arg_dict.get("flags", [])
            if not flags:
                return None

            # Get the main flag (usually the long form)
            main_flag = next(
                (flag for flag in flags if flag.startswith("--")), flags[0]
            )
            name = main_flag.lstrip("-").replace("-", "_")

            # Extract other properties
            description = arg_dict.get("help", "")
            default_value = arg_dict.get("default")
            choices = arg_dict.get("choices")
            required = arg_dict.get("required", False)

            # Infer type from argument configuration
            param_type = self._infer_arg_type(arg_dict)

            return ConfigOption(
                name=name,
                description=description,
                type=param_type,
                default=default_value,
                choices=choices,
                required=required,
                section=section,
            )

        except Exception as e:
            logger.warning(f"Failed to create config option from arg: {e}")
            return None

    def _infer_arg_type(self, arg_dict: dict) -> str:
        """Infer the type of an argument from its configuration."""
        action = arg_dict.get("action", "")

        if action in ["store_true", "store_false"]:
            return "bool"
        elif action == "append":
            return "list"
        elif "type" in arg_dict:
            arg_type = arg_dict["type"]
            if hasattr(arg_type, "__name__"):
                return arg_type.__name__
            else:
                return str(arg_type)
        elif "choices" in arg_dict:
            return "choice"
        elif arg_dict.get("default") is not None:
            return type(arg_dict["default"]).__name__
        else:
            return "str"


class OptionsDiscovery:
    """Main class for discovering all immunopipe configuration options."""

    def __init__(self):
        self.pipeline_discovery = PipelineOptionsDiscovery()
        self.process_discovery = ProcessDiscovery()
        self.gbatch_discovery = GbatchOptionsDiscovery()

    def get_all_options(self) -> Dict[str, Any]:
        """Get all available options categorized by type."""
        return {
            "pipeline": self.pipeline_discovery.get_pipeline_options(),
            "processes": self.process_discovery.get_processes(),
            "gbatch": self.gbatch_discovery.get_gbatch_options(),
        }

    def get_pipeline_options(self) -> Dict[str, ConfigOption]:
        """Get pipeline options."""
        return self.pipeline_discovery.get_pipeline_options()

    def get_process_options(self) -> Dict[str, Dict[str, Any]]:
        """Get process options."""
        return self.process_discovery.get_processes()

    def get_gbatch_options(self) -> Dict[str, ConfigOption]:
        """Get gbatch options."""
        return self.gbatch_discovery.get_gbatch_options()


_CONFIG_SURFACE: Optional[Dict[str, Any]] = None


def get_config_surface() -> Dict[str, Any]:
    """Get the configuration surface of the installed pipeline.

    The result is cached, as building the pipeline is expensive.

    Returns:
        A mapping with:
        - `pipeline_options`: the names of the pipeline-level options
        - `processes`: a mapping of process name to the names of its `envs`,
          `inputs` and `outputs`
    """
    global _CONFIG_SURFACE

    if _CONFIG_SURFACE is None:
        pipeline_options = set(PipelineOptionsDiscovery().get_pipeline_options())
        # Some pipeline options are not exposed by `PIPEN_ARGS`
        pipeline_options.update(PIPEN_CONFIG)
        _CONFIG_SURFACE = {
            "pipeline_options": pipeline_options,
            "processes": ProcessDiscovery().get_config_surface(),
        }

    return _CONFIG_SURFACE
