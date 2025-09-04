"""MCP tools for immunopipe configuration generation."""

from typing import Dict, Any, List, Optional
import logging
from dataclasses import dataclass

from .options import OptionsDiscovery
from .toml_generator import TOMLGenerator, ConfigTemplateGenerator
from .doc_extractor import ProcessDocumentationExtractor

logger = logging.getLogger(__name__)


@dataclass
class MCPToolResult:
    """Result from an MCP tool execution."""

    success: bool
    content: Any
    message: str = ""
    tool_calls: Optional[List[Dict[str, Any]]] = None


class ImmunopipeConfigTools:
    """Hierarchical MCP tools for immunopipe configuration generation."""

    def __init__(self):
        self.options_discovery = OptionsDiscovery()
        self.toml_generator = TOMLGenerator()
        self.template_generator = ConfigTemplateGenerator(self.toml_generator)
        self.doc_extractor = ProcessDocumentationExtractor()

        # Tool registry for hierarchical calls
        self.tools = {
            # Discovery tools
            "list_pipeline_options": self.list_pipeline_options,
            "list_processes": self.list_processes,
            "list_gbatch_options": self.list_gbatch_options,
            "get_process_details": self.get_process_details,
            # Generation tools
            "generate_pipeline_config": self.generate_pipeline_config,
            "generate_process_config": self.generate_process_config,
            "generate_gbatch_config": self.generate_gbatch_config,
            "generate_full_config": self.generate_full_config,
            # Template tools
            "generate_basic_template": self.generate_basic_template,
            "generate_tcr_template": self.generate_tcr_template,
            "generate_gbatch_template": self.generate_gbatch_template,
            # Manipulation tools
            "merge_configs": self.merge_configs,
            "validate_config": self.validate_config,
            "format_config": self.format_config,
            # Assistant tools
            "configure_immunopipe": self.configure_immunopipe,
            "help_generate_config": self.help_generate_config,
            "suggest_processes": self.suggest_processes,
        }

    def get_available_tools(self) -> List[Dict[str, Any]]:
        """Get list of available MCP tools."""
        return [
            {
                "name": "list_pipeline_options",
                "description": "List all available pipeline configuration options",
                "inputSchema": {"type": "object", "properties": {}, "required": []},
            },
            {
                "name": "list_processes",
                "description": "List all available immunopipe processes",
                "inputSchema": {"type": "object", "properties": {}, "required": []},
            },
            {
                "name": "list_gbatch_options",
                "description": "List all Google Batch configuration options",
                "inputSchema": {"type": "object", "properties": {}, "required": []},
            },
            {
                "name": "get_process_details",
                "description": "Get detailed information about a specific process",
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "process_name": {
                            "type": "string",
                            "description": "Name of the process to get details for",
                        }
                    },
                    "required": ["process_name"],
                },
            },
            {
                "name": "generate_pipeline_config",
                "description": "Generate pipeline-level configuration",
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "options": {
                            "type": "object",
                            "description": "Pipeline options to include",
                        }
                    },
                    "required": ["options"],
                },
            },
            {
                "name": "generate_process_config",
                "description": "Generate configuration for a specific process",
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "process_name": {
                            "type": "string",
                            "description": "Name of the process",
                        },
                        "config": {
                            "type": "object",
                            "description": "Process configuration options",
                        },
                    },
                    "required": ["process_name", "config"],
                },
            },
            {
                "name": "generate_gbatch_config",
                "description": "Generate Google Batch configuration",
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "options": {
                            "type": "object",
                            "description": "Google Batch options",
                        }
                    },
                    "required": ["options"],
                },
            },
            {
                "name": "generate_full_config",
                "description": (
                    "Generate a complete immunopipe TOML configuration file "
                    "for single-cell RNA-seq and TCR/BCR analysis. "
                    "Use this tool when users ask to create immunopipe configurations, "
                    "set clustering parameters, configure analysis pipelines, "
                    "or generate TOML config files."
                ),
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "pipeline_options": {
                            "type": "object",
                            "description": (
                                "Pipeline-level options " "(name, outdir, forks, etc.)"
                            ),
                        },
                        "processes": {
                            "type": "object",
                            "description": (
                                "Process configurations (SeuratClustering, "
                                "TCRClustering, etc.)"
                            ),
                        },
                        "gbatch_options": {
                            "type": "object",
                            "description": "Google Batch options for cloud execution",
                        },
                        "description": {
                            "type": "string",
                            "description": "Description for the configuration file",
                        },
                    },
                    "required": [],
                },
            },
            {
                "name": "generate_basic_template",
                "description": "Generate a basic configuration template",
                "inputSchema": {"type": "object", "properties": {}, "required": []},
            },
            {
                "name": "generate_tcr_template",
                "description": "Generate a template for TCR analysis",
                "inputSchema": {"type": "object", "properties": {}, "required": []},
            },
            {
                "name": "generate_gbatch_template",
                "description": "Generate a template for Google Batch execution",
                "inputSchema": {"type": "object", "properties": {}, "required": []},
            },
            {
                "name": "merge_configs",
                "description": "Merge two configuration files or strings",
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "base_config": {
                            "type": "string",
                            "description": "Base configuration (TOML string)",
                        },
                        "new_config": {
                            "type": "string",
                            "description": "New configuration to merge (TOML string)",
                        },
                    },
                    "required": ["base_config", "new_config"],
                },
            },
            {
                "name": "validate_config",
                "description": "Validate a TOML configuration",
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "config_content": {
                            "type": "string",
                            "description": "TOML configuration content to validate",
                        }
                    },
                    "required": ["config_content"],
                },
            },
            {
                "name": "format_config",
                "description": "Format and enhance a TOML configuration",
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "config_content": {
                            "type": "string",
                            "description": "TOML configuration content to format",
                        },
                        "add_comments": {
                            "type": "boolean",
                            "description": "Whether to add helpful comments",
                        },
                    },
                    "required": ["config_content"],
                },
            },
            {
                "name": "configure_immunopipe",
                "description": (
                    "Generate immunopipe configuration based on natural "
                    "language requirements. Automatically identifies required "
                    "processes, extracts their documentation, and generates "
                    "appropriate TOML configuration. Use this for any immunopipe "
                    "configuration request (clustering, TCR analysis, process "
                    "parameters, etc.)."
                ),
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "requirements": {
                            "type": "string",
                            "description": (
                                "Natural language description of what you want to "
                                "configure (e.g., 'clustering with resolution 0.5', "
                                "'TCR analysis with clustering', 'set SeuratPreparing "
                                "QC parameters')"
                            ),
                        },
                        "specific_processes": {
                            "type": "array",
                            "items": {"type": "string"},
                            "description": (
                                "Optional: specific process names to configure "
                                "(if not provided, will be inferred from requirements)"
                            ),
                        },
                        "parameters": {
                            "type": "object",
                            "description": "Optional: specific parameter values to set",
                        },
                    },
                    "required": ["requirements"],
                },
            },
            {
                "name": "help_generate_config",
                "description": (
                    "Get help and suggestions for generating configuration based on "
                    "natural language description. Use this when users ask for help "
                    "with immunopipe configuration or describe what they "
                    "want to analyze."
                ),
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "description": {
                            "type": "string",
                            "description": (
                                "Natural language description of what you want to "
                                "configure (e.g., 'I want to analyze TCR data', "
                                "'I need clustering analysis')"
                            ),
                        }
                    },
                    "required": ["description"],
                },
            },
            {
                "name": "suggest_processes",
                "description": (
                    "Suggest processes based on analysis type or requirements"
                ),
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "analysis_type": {
                            "type": "string",
                            "description": (
                                "Type of analysis (e.g., 'tcr', 'bcr', 'clustering', "
                                "'differential')"
                            ),
                        }
                    },
                    "required": ["analysis_type"],
                },
            },
        ]

    async def execute_tool(
        self, tool_name: str, parameters: Dict[str, Any]
    ) -> MCPToolResult:
        """Execute a specific MCP tool."""
        try:
            if tool_name not in self.tools:
                return MCPToolResult(
                    success=False, content=None, message=f"Unknown tool: {tool_name}"
                )

            tool_func = self.tools[tool_name]
            result = await tool_func(**parameters)

            if isinstance(result, MCPToolResult):
                return result
            else:
                return MCPToolResult(success=True, content=result)

        except Exception as e:
            logger.error(f"Error executing tool {tool_name}: {e}")
            return MCPToolResult(
                success=False, content=None, message=f"Tool execution error: {str(e)}"
            )

    # Discovery tools
    async def list_pipeline_options(self) -> MCPToolResult:
        """List all available pipeline options."""
        try:
            options = self.options_discovery.get_pipeline_options()
            formatted_options = {}

            for name, option in options.items():
                formatted_options[name] = {
                    "description": option.description,
                    "type": option.type,
                    "default": option.default,
                    "required": option.required,
                }
                if option.choices:
                    formatted_options[name]["choices"] = option.choices

            return MCPToolResult(
                success=True,
                content=formatted_options,
                message="Retrieved pipeline options successfully",
            )
        except Exception as e:
            return MCPToolResult(
                success=False,
                content=None,
                message=f"Failed to retrieve pipeline options: {e}",
            )

    async def list_processes(self) -> MCPToolResult:
        """List all available processes."""
        try:
            processes = self.options_discovery.get_process_options()
            formatted_processes = {}

            for name, process_info in processes.items():
                formatted_processes[name] = {
                    "description": process_info["description"],
                    "envs": process_info.get("envs", {}),
                    "requires": process_info.get("requires", []),
                }

            return MCPToolResult(
                success=True,
                content=formatted_processes,
                message="Retrieved processes successfully",
            )
        except Exception as e:
            return MCPToolResult(
                success=False,
                content=None,
                message=f"Failed to retrieve processes: {e}",
            )

    async def list_gbatch_options(self) -> MCPToolResult:
        """List all Google Batch options."""
        try:
            options = self.options_discovery.get_gbatch_options()
            formatted_options = {}

            for name, option in options.items():
                formatted_options[name] = {
                    "description": option.description,
                    "type": option.type,
                    "default": option.default,
                }

            return MCPToolResult(
                success=True,
                content=formatted_options,
                message="Retrieved Google Batch options successfully",
            )
        except Exception as e:
            return MCPToolResult(
                success=False,
                content=None,
                message=f"Failed to retrieve Google Batch options: {e}",
            )

    async def get_process_details(self, process_name: str) -> MCPToolResult:
        """Get detailed information about a specific process."""
        try:
            processes = self.options_discovery.get_process_options()

            if process_name not in processes:
                return MCPToolResult(
                    success=False,
                    content=None,
                    message=f"Process '{process_name}' not found",
                )

            process_info = processes[process_name]

            return MCPToolResult(
                success=True,
                content=process_info,
                message=f"Retrieved details for process '{process_name}'",
            )
        except Exception as e:
            return MCPToolResult(
                success=False,
                content=None,
                message=f"Failed to get process details: {e}",
            )

    # Generation tools
    async def generate_pipeline_config(self, options: Dict[str, Any]) -> MCPToolResult:
        """Generate pipeline configuration."""
        try:
            config_content = self.toml_generator.generate_pipeline_section(options)

            return MCPToolResult(
                success=True,
                content=config_content,
                message="Generated pipeline configuration successfully",
            )
        except Exception as e:
            return MCPToolResult(
                success=False,
                content=None,
                message=f"Failed to generate pipeline config: {e}",
            )

    async def generate_process_config(
        self, process_name: str, config: Dict[str, Any]
    ) -> MCPToolResult:
        """Generate process configuration."""
        try:
            config_content = self.toml_generator.generate_process_section(
                process_name, config
            )

            return MCPToolResult(
                success=True,
                content=config_content,
                message=f"Generated configuration for process '{process_name}'",
            )
        except Exception as e:
            return MCPToolResult(
                success=False,
                content=None,
                message=f"Failed to generate process config: {e}",
            )

    async def generate_gbatch_config(self, options: Dict[str, Any]) -> MCPToolResult:
        """Generate Google Batch configuration."""
        try:
            config_content = self.toml_generator.generate_gbatch_section(options)

            return MCPToolResult(
                success=True,
                content=config_content,
                message="Generated Google Batch configuration successfully",
            )
        except Exception as e:
            return MCPToolResult(
                success=False,
                content=None,
                message=f"Failed to generate gbatch config: {e}",
            )

    async def generate_full_config(
        self,
        pipeline_options: Dict[str, Any] = None,
        processes: Dict[str, Dict[str, Any]] = None,
        gbatch_options: Dict[str, Any] = None,
        description: str = "Immunopipe configuration",
    ) -> MCPToolResult:
        """Generate a complete configuration file."""
        try:
            config_content = self.toml_generator.generate_full_config(
                pipeline_options=pipeline_options,
                process_configs=processes,
                gbatch_options=gbatch_options,
                description=description,
            )

            return MCPToolResult(
                success=True,
                content=config_content,
                message="Generated complete configuration successfully",
            )
        except Exception as e:
            return MCPToolResult(
                success=False,
                content=None,
                message=f"Failed to generate full config: {e}",
            )

    # Template tools
    async def generate_basic_template(self) -> MCPToolResult:
        """Generate a basic configuration template."""
        try:
            template = self.template_generator.generate_basic_template()

            return MCPToolResult(
                success=True,
                content=template,
                message="Generated basic template successfully",
            )
        except Exception as e:
            return MCPToolResult(
                success=False,
                content=None,
                message=f"Failed to generate basic template: {e}",
            )

    async def generate_tcr_template(self) -> MCPToolResult:
        """Generate a TCR analysis template."""
        try:
            template = self.template_generator.generate_tcr_analysis_template()

            return MCPToolResult(
                success=True,
                content=template,
                message="Generated TCR analysis template successfully",
            )
        except Exception as e:
            return MCPToolResult(
                success=False,
                content=None,
                message=f"Failed to generate TCR template: {e}",
            )

    async def generate_gbatch_template(self) -> MCPToolResult:
        """Generate a Google Batch template."""
        try:
            template = self.template_generator.generate_gbatch_template()

            return MCPToolResult(
                success=True,
                content=template,
                message="Generated Google Batch template successfully",
            )
        except Exception as e:
            return MCPToolResult(
                success=False,
                content=None,
                message=f"Failed to generate gbatch template: {e}",
            )

    # Manipulation tools
    async def merge_configs(self, base_config: str, new_config: str) -> MCPToolResult:
        """Merge two configurations."""
        try:
            merged_config = self.toml_generator.merge_configs(base_config, new_config)

            return MCPToolResult(
                success=True,
                content=merged_config,
                message="Merged configurations successfully",
            )
        except Exception as e:
            return MCPToolResult(
                success=False, content=None, message=f"Failed to merge configs: {e}"
            )

    async def validate_config(self, config_content: str) -> MCPToolResult:
        """Validate configuration."""
        try:
            is_valid, errors = self.toml_generator.validate_config(config_content)

            return MCPToolResult(
                success=True,
                content={"valid": is_valid, "errors": errors},
                message="Validation completed",
            )
        except Exception as e:
            return MCPToolResult(
                success=False, content=None, message=f"Failed to validate config: {e}"
            )

    async def format_config(
        self, config_content: str, add_comments: bool = True
    ) -> MCPToolResult:
        """Format configuration."""
        try:
            formatted_config = self.toml_generator.format_config(
                config_content, add_comments
            )

            return MCPToolResult(
                success=True,
                content=formatted_config,
                message="Formatted configuration successfully",
            )
        except Exception as e:
            return MCPToolResult(
                success=False, content=None, message=f"Failed to format config: {e}"
            )

    # Assistant tools
    async def configure_immunopipe(
        self,
        requirements: str,
        specific_processes: List[str] = None,
        parameters: Dict[str, Any] = None,
    ) -> MCPToolResult:
        """Generate immunopipe configuration based on natural language requirements."""
        try:
            if parameters is None:
                parameters = {}

            # Get all available processes
            all_processes = self.options_discovery.get_process_options()

            # Analyze requirements to identify needed processes
            required_processes = self._analyze_requirements(
                requirements, all_processes, specific_processes
            )

            # Generate configuration for identified processes
            process_configs = {}
            pipeline_options = {"name": "immunopipe_analysis", "outdir": "./output"}

            for process_name in required_processes:
                process_config = self._generate_process_config(
                    process_name, requirements, parameters, all_processes
                )
                # Include process even if no specific config (enables it in pipeline)
                process_configs[process_name] = process_config if process_config else {}

            # Generate the complete configuration
            config_content = self.toml_generator.generate_full_config(
                pipeline_options=pipeline_options,
                process_configs=process_configs,
                description=f"Immunopipe configuration for: {requirements}",
            )

            return MCPToolResult(
                success=True,
                content=config_content,
                message=(
                    "Generated immunopipe configuration for processes: "
                    f"{', '.join(required_processes)}"
                ),
            )

        except Exception as e:
            return MCPToolResult(
                success=False,
                content=None,
                message=f"Failed to configure immunopipe: {e}",
            )

    def _analyze_requirements(
        self,
        requirements: str,
        all_processes: Dict,
        specific_processes: List[str] = None,
    ) -> List[str]:
        """Analyze natural language requirements to identify needed processes."""
        if specific_processes:
            return specific_processes

        required = []
        requirements_lower = requirements.lower()

        # Always include SampleInfo as it's the entry point
        required.append("SampleInfo")

        # Determine if TCR/BCR data is involved (affects process selection)
        has_vdj_data = any(
            keyword in requirements_lower
            for keyword in [
                "tcr",
                "tcell",
                "t cell",
                "t-cell",
                "bcr",
                "bcell",
                "b cell",
                "b-cell",
                "vdj",
            ]
        )

        # Analysis type detection
        if any(
            keyword in requirements_lower
            for keyword in ["tcr", "tcell", "t cell", "t-cell"]
        ):
            required.extend(["ScRepLoading", "SeuratPreparing"])
            # Use clustering of all cells first, then T cell selection
            if "SeuratClusteringOfAllCells" in all_processes:
                required.append("SeuratClusteringOfAllCells")
            required.extend(
                [
                    "TOrBCellSelection",
                    "SeuratClustering",
                    "TCRClustering",
                    "ClonalStats",
                ]
            )
        elif any(
            keyword in requirements_lower
            for keyword in ["bcr", "bcell", "b cell", "b-cell"]
        ):
            required.extend(["ScRepLoading", "SeuratPreparing"])
            # Use clustering of all cells first, then B cell selection
            if "SeuratClusteringOfAllCells" in all_processes:
                required.append("SeuratClusteringOfAllCells")
            required.extend(["TOrBCellSelection", "SeuratClustering", "ClonalStats"])
        elif any(
            keyword in requirements_lower for keyword in ["cluster", "clustering"]
        ):
            required.append("SeuratPreparing")
            # Choose appropriate clustering based on context
            clustering_process = self._choose_clustering_process(
                requirements_lower, all_processes, has_vdj_data
            )
            required.append(clustering_process)
        elif any(
            keyword in requirements_lower
            for keyword in ["marker", "differential", "deg"]
        ):
            required.append("SeuratPreparing")
            clustering_process = self._choose_clustering_process(
                requirements_lower, all_processes, has_vdj_data
            )
            required.extend([clustering_process, "MarkersFinder"])
        else:
            # Basic RNA-seq analysis
            required.append("SeuratPreparing")
            clustering_process = self._choose_clustering_process(
                requirements_lower, all_processes, has_vdj_data
            )
            required.append(clustering_process)

        # Process-specific detection for explicitly mentioned processes
        for process_name in all_processes.keys():
            if process_name.lower() in requirements_lower:
                if process_name not in required:
                    required.append(process_name)

        # Resolve conflicts and choose optimal processes
        required = self._resolve_process_conflicts(
            required, requirements_lower, all_processes, has_vdj_data
        )

        # Filter to only include processes that exist
        return [p for p in required if p in all_processes]

    def _choose_clustering_process(
        self, requirements_lower: str, all_processes: Dict, has_vdj_data: bool
    ) -> str:
        """Choose the appropriate clustering process based on context."""
        # If TCR/BCR data is involved and we have the "OfAllCells" variant, use it
        if has_vdj_data and "SeuratClusteringOfAllCells" in all_processes:
            return "SeuratClusteringOfAllCells"

        # If explicitly asking for all cells clustering
        if "all cells" in requirements_lower or "allcells" in requirements_lower:
            return (
                "SeuratClusteringOfAllCells"
                if "SeuratClusteringOfAllCells" in all_processes
                else "SeuratClustering"
            )

        # Default to standard clustering
        return "SeuratClustering"

    def _resolve_process_conflicts(
        self,
        required: List[str],
        requirements_lower: str,
        all_processes: Dict,
        has_vdj_data: bool,
    ) -> List[str]:
        """Resolve conflicts between similar processes and choose the best ones."""
        resolved = []

        for process in required:
            # Handle clustering process conflicts
            if (
                process == "SeuratClustering"
                and "SeuratClusteringOfAllCells" in required
            ):
                # Keep both - OfAllCells comes first,
                # then regular clustering after cell selection
                if process not in resolved:
                    resolved.append(process)
            elif process == "SeuratClusteringOfAllCells":
                if process not in resolved:
                    resolved.append(process)
            # Handle marker finder conflicts
            elif process == "MarkersFinder" and "ClusterMarkersOfAllCells" in required:
                # Choose based on context -
                # if we have all cells clustering, use all cells markers first
                if "SeuratClusteringOfAllCells" in required:
                    if "ClusterMarkersOfAllCells" not in resolved:
                        resolved.append("ClusterMarkersOfAllCells")
                    if process not in resolved:
                        resolved.append(process)
                else:
                    if process not in resolved:
                        resolved.append(process)
            elif process == "ClusterMarkersOfAllCells":
                if process not in resolved:
                    resolved.append(process)
            # Handle top expressing genes conflicts
            elif (
                process == "TopExpressingGenes"
                and "TopExpressingGenesOfAllCells" in required
            ):
                if "SeuratClusteringOfAllCells" in required:
                    if "TopExpressingGenesOfAllCells" not in resolved:
                        resolved.append("TopExpressingGenesOfAllCells")
                    if process not in resolved:
                        resolved.append(process)
                else:
                    if process not in resolved:
                        resolved.append(process)
            elif process == "TopExpressingGenesOfAllCells":
                if process not in resolved:
                    resolved.append(process)
            else:
                if process not in resolved:
                    resolved.append(process)

        return resolved

    def _generate_process_config(
        self,
        process_name: str,
        requirements: str,
        parameters: Dict,
        all_processes: Dict,
    ) -> Dict[str, Any]:
        """Generate configuration for a specific process based on requirements."""
        try:
            # Get the process class to extract documentation
            process_class = self._get_process_class(process_name)
            if not process_class:
                logger.warning(f"Could not find process class for {process_name}")
                return self._fallback_process_config(
                    process_name, requirements, parameters
                )

            # Use generic documentation extractor
            extracted_config = self.doc_extractor.extract_parameters_from_requirements(
                requirements, process_class
            )

            # Add any user-specified parameters for this process
            if parameters.get(process_name):
                if "envs" not in extracted_config:
                    extracted_config["envs"] = {}
                extracted_config["envs"].update(parameters[process_name])

            return extracted_config if extracted_config else {}

        except Exception as e:
            logger.warning(f"Error extracting config for {process_name}: {e}")
            return self._fallback_process_config(process_name, requirements, parameters)

    def _get_process_class(self, process_name: str):
        """Get the process class by name."""
        try:
            # Try to import the process class dynamically
            import immunopipe.processes

            if hasattr(immunopipe.processes, process_name):
                return getattr(immunopipe.processes, process_name)

            # Some processes might be conditional, try to get from pipeline
            from immunopipe.pipeline import Immunopipe

            pipeline = Immunopipe()

            for proc in pipeline.procs:
                if proc.name == process_name:
                    return proc.__class__

            return None

        except Exception as e:
            logger.error(f"Error getting process class for {process_name}: {e}")
            return None

    def _fallback_process_config(
        self, process_name: str, requirements: str, parameters: Dict
    ) -> Dict[str, Any]:
        """Fallback configuration generation for when documentation extraction fails."""
        config = {}
        requirements_lower = requirements.lower()

        # Basic fallback for common processes
        if process_name in ["SeuratClustering", "SeuratClusteringOfAllCells"]:
            envs = {}

            # Extract resolution parameter
            import re

            resolution_match = re.search(r"resolution.*?([0-9.]+)", requirements_lower)
            if resolution_match:
                envs["FindClusters"] = {"resolution": float(resolution_match.group(1))}
            elif "cluster" in requirements_lower:
                envs["FindClusters"] = {"resolution": 0.8}  # Default

            if envs:
                config["envs"] = envs

        elif process_name == "TOrBCellSelection":
            envs = {}
            if "tcr" in requirements_lower or "t cell" in requirements_lower:
                envs["cell_type"] = "T"
            elif "bcr" in requirements_lower or "b cell" in requirements_lower:
                envs["cell_type"] = "B"
            if envs:
                config["envs"] = envs

        # Add user-specified parameters
        if parameters.get(process_name):
            if "envs" not in config:
                config["envs"] = {}
            config["envs"].update(parameters[process_name])

        return config if config else {}

    async def help_generate_config(self, description: str) -> MCPToolResult:
        """Generate configuration based on natural language description."""
        try:
            # This would typically call other tools based on the description
            suggestions = []
            tools_to_call = []

            # Simple keyword-based suggestions
            description_lower = description.lower()

            if any(
                keyword in description_lower for keyword in ["tcr", "tcell", "t cell"]
            ):
                suggestions.append(
                    "Consider using TCR analysis processes: TOrBCellSelection, "
                    "TCRClustering, ClonalStats"
                )
                tools_to_call.append(
                    {
                        "tool": "suggest_processes",
                        "parameters": {"analysis_type": "tcr"},
                    }
                )

            if any(
                keyword in description_lower for keyword in ["bcr", "bcell", "b cell"]
            ):
                suggestions.append("Consider using BCR analysis processes")
                tools_to_call.append(
                    {
                        "tool": "suggest_processes",
                        "parameters": {"analysis_type": "bcr"},
                    }
                )

            if any(
                keyword in description_lower for keyword in ["cluster", "clustering"]
            ):
                suggestions.append(
                    "Consider using clustering processes: SeuratClustering, "
                    "SeuratSubClustering"
                )
                tools_to_call.append(
                    {
                        "tool": "suggest_processes",
                        "parameters": {"analysis_type": "clustering"},
                    }
                )

            if any(
                keyword in description_lower for keyword in ["google", "batch", "cloud"]
            ):
                suggestions.append("Consider using Google Batch for cloud execution")
                tools_to_call.append(
                    {"tool": "generate_gbatch_template", "parameters": {}}
                )

            if not suggestions:
                suggestions.append("For basic analysis, start with the basic template")
                tools_to_call.append(
                    {"tool": "generate_basic_template", "parameters": {}}
                )

            return MCPToolResult(
                success=True,
                content={
                    "suggestions": suggestions,
                    "recommended_tools": tools_to_call,
                },
                message="Generated configuration suggestions",
                tool_calls=tools_to_call,
            )

        except Exception as e:
            return MCPToolResult(
                success=False, content=None, message=f"Failed to generate help: {e}"
            )

    async def suggest_processes(self, analysis_type: str) -> MCPToolResult:
        """Suggest processes based on analysis type."""
        try:
            processes_map = {
                "tcr": [
                    "SampleInfo",
                    "ScRepLoading",
                    "SeuratPreparing",
                    "TOrBCellSelection",
                    "SeuratClustering",
                    "TCRClustering",
                    "ClonalStats",
                    "TESSA",
                ],
                "bcr": [
                    "SampleInfo",
                    "ScRepLoading",
                    "SeuratPreparing",
                    "TOrBCellSelection",
                    "SeuratClustering",
                    "ClonalStats",
                ],
                "clustering": [
                    "SampleInfo",
                    "SeuratPreparing",
                    "SeuratClustering",
                    "SeuratSubClustering",
                    "CellTypeAnnotation",
                    "ClusterMarkers",
                ],
                "differential": [
                    "SampleInfo",
                    "SeuratPreparing",
                    "SeuratClustering",
                    "MarkersFinder",
                    "PseudoBulkDEG",
                ],
                "basic": [
                    "SampleInfo",
                    "SeuratPreparing",
                    "SeuratClustering",
                    "SeuratClusterStats",
                ],
            }

            suggested_processes = processes_map.get(
                analysis_type.lower(), processes_map["basic"]
            )

            return MCPToolResult(
                success=True,
                content={
                    "analysis_type": analysis_type,
                    "suggested_processes": suggested_processes,
                    "description": (
                        f"Recommended processes for {analysis_type} analysis"
                    ),
                },
                message=f"Generated process suggestions for {analysis_type} analysis",
            )

        except Exception as e:
            return MCPToolResult(
                success=False, content=None, message=f"Failed to suggest processes: {e}"
            )
