"""Generic documentation-based parameter extractor using pipen-annotate."""

from typing import Dict, Any, Optional
import re
import logging
from pipen_annotate import annotate

logger = logging.getLogger(__name__)


class ProcessDocumentationExtractor:
    """Extract configuration parameters from process docs."""

    def __init__(self):
        self._process_annotations_cache: Dict[str, Any] = {}

    def get_process_annotations(self, process_class) -> Dict[str, Any]:
        """Get annotated documentation for a process class."""
        process_name = process_class.__name__

        if process_name in self._process_annotations_cache:
            return self._process_annotations_cache[process_name]

        try:
            annotations = annotate(process_class)
            self._process_annotations_cache[process_name] = dict(annotations)
            return self._process_annotations_cache[process_name]
        except Exception as e:
            logger.error(f"Failed to annotate process {process_name}: {e}")
            return {}

    def extract_env_parameters(self, process_class) -> Dict[str, Any]:
        """Extract environment parameters from process documentation."""
        annotations = self.get_process_annotations(process_class)

        # pipen-annotate returns capitalized keys
        env_key = 'Envs' if 'Envs' in annotations else 'envs'
        if env_key not in annotations:
            return {}

        return dict(annotations[env_key])

    def extract_parameters_from_requirements(
        self,
        requirements: str,
        process_class
    ) -> Dict[str, Any]:
        """Extract configuration parameters based on requirements text."""
        env_params = self.extract_env_parameters(process_class)

        if not env_params:
            return {}

        extracted_config = {}
        requirements_lower = requirements.lower()

        # Process each environment parameter
        for param_name, param_info in env_params.items():
            param_config = self._extract_parameter_config(
                param_name, param_info, requirements_lower
            )

            if param_config:
                extracted_config.update(param_config)

        return extracted_config

    def _extract_parameter_config(
        self,
        param_name: str,
        param_info: Dict[str, Any],
        requirements_lower: str
    ) -> Dict[str, Any]:
        """Extract configuration for a specific parameter."""
        config = {}

        # Handle namespace parameters (like RunUMAP, FindClusters)
        # In pipen-annotate, namespace is indicated by attrs.ns = true
        attrs = param_info.get('attrs', {})
        if attrs.get('ns'):
            ns_config = self._extract_namespace_config(
                param_name, param_info, requirements_lower
            )
            if ns_config:
                config[param_name] = ns_config

        # Handle regular parameters
        else:
            param_value = self._extract_parameter_value(
                param_name, param_info, requirements_lower
            )
            if param_value is not None:
                config[param_name] = param_value

        return config

    def _extract_namespace_config(
        self,
        ns_name: str,
        ns_info: Dict[str, Any],
        requirements_lower: str
    ) -> Dict[str, Any]:
        """Extract configuration for namespace parameters."""
        ns_config = {}

        # Parse nested parameters from the help text
        help_text = ns_info.get('help', '')
        nested_params = self._parse_nested_params_from_help(help_text)

        for nested_name, nested_info in nested_params.items():
            nested_value = self._extract_parameter_value(
                nested_name, nested_info, requirements_lower, parent_ns=ns_name
            )
            if nested_value is not None:
                ns_config[nested_name] = nested_value

        return ns_config if ns_config else None

    def _parse_nested_params_from_help(
        self, help_text: str
    ) -> Dict[str, Dict[str, Any]]:
        """Parse nested parameters from help text."""
        nested_params = {}

        # Look for patterns like "dims (type=int): Number of dimensions"
        param_pattern = r'(\w+)\s*\(type=(\w+)\):\s*([^.]+)'

        for match in re.finditer(param_pattern, help_text):
            param_name = match.group(1)
            param_type = match.group(2)
            description = match.group(3).strip()

            nested_params[param_name] = {
                'type': param_type,
                'description': description
            }

        # Also look for parameters without explicit type
        simple_pattern = r'(\w+):\s*([^.]+(?:\.|$))'

        for match in re.finditer(simple_pattern, help_text):
            param_name = match.group(1)
            description = match.group(2).strip()

            # Skip if already found with type info
            if param_name not in nested_params:
                nested_params[param_name] = {
                    'type': 'str',
                    'description': description
                }

        return nested_params

    def _extract_parameter_value(
        self,
        param_name: str,
        param_info: Dict[str, Any],
        requirements_lower: str,
        parent_ns: Optional[str] = None
    ) -> Any:
        """Extract a specific parameter value from requirements text."""
        param_name_lower = param_name.lower()

        # Get parameter type from attrs or direct type info
        if isinstance(param_info, dict):
            if 'attrs' in param_info:
                param_type = param_info['attrs'].get('type', 'str')
            else:
                param_type = param_info.get('type', 'str')
        else:
            param_type = 'str'

        # Look for explicit parameter mentions
        patterns = [
            rf'{re.escape(param_name_lower)}\s*[=:]\s*([^\s,;]+)',
            rf'{re.escape(param_name_lower)}\s+([0-9.]+)',
            rf'set\s+{re.escape(param_name_lower)}\s+to\s+([^\s,;]+)',
        ]

        for pattern in patterns:
            match = re.search(pattern, requirements_lower)
            if match:
                value_str = match.group(1)
                return self._convert_value(value_str, param_type)

        # Apply intelligent suggestions
        return self._suggest_parameter_value(
            param_name_lower, param_type, requirements_lower, parent_ns
        )

    def _suggest_parameter_value(
        self,
        param_name: str,
        param_type: str,
        requirements_lower: str,
        parent_ns: Optional[str] = None
    ) -> Any:
        """Suggest parameter value based on context."""
        # Resolution parameter handling
        if param_name == 'resolution':
            resolution_match = re.search(
                r'resolution.*?([0-9.]+)', requirements_lower
            )
            if resolution_match:
                return float(resolution_match.group(1))
            elif 'high resolution' in requirements_lower:
                return 1.2
            elif 'low resolution' in requirements_lower:
                return 0.3
            elif 'fine' in requirements_lower:
                return 1.0
            elif 'cluster' in requirements_lower:
                return 0.8

        # Algorithm parameter handling
        elif param_name == 'algorithm':
            if 'leiden' in requirements_lower:
                return 'Leiden'
            elif 'louvain' in requirements_lower:
                return 'Louvain'

        # Dimensions parameter handling
        elif param_name == 'dims':
            dims_match = re.search(
                r'(\d+)\s*(?:dimensions?|dims?|pcs?)', requirements_lower
            )
            if dims_match:
                return int(dims_match.group(1))
            elif parent_ns == 'RunUMAP':
                return 30

        # QC parameters
        elif param_name in ['min_cells', 'min_features']:
            if 'qc' in requirements_lower:
                if param_name == 'min_cells':
                    return 3
                elif param_name == 'min_features':
                    return 200

        # Cell type selection
        elif param_name == 'cell_type':
            if 'tcr' in requirements_lower or 't cell' in requirements_lower:
                return 'T'
            elif 'bcr' in requirements_lower or 'b cell' in requirements_lower:
                return 'B'

        # Mitochondrial pattern
        elif param_name == 'mito_pattern':
            if 'mitochondrial' in requirements_lower:
                return '^MT-'

        # QC filters
        elif param_name == 'cell_qc':
            if 'qc' in requirements_lower:
                return "nFeature_RNA > 200 & percent.mt < 5"

        return None

    def _convert_value(self, value_str: str, param_type: str) -> Any:
        """Convert string value to appropriate type."""
        value_str = value_str.strip()

        if param_type == 'int':
            try:
                return int(float(value_str))
            except ValueError:
                return value_str
        elif param_type == 'float':
            try:
                return float(value_str)
            except ValueError:
                return value_str
        elif param_type == 'bool':
            if value_str.lower() in ('true', 'yes', '1'):
                return True
            elif value_str.lower() in ('false', 'no', '0'):
                return False
            return value_str
        else:
            return value_str

    def get_all_process_parameters(self, process_class) -> Dict[str, Any]:
        """Get all available parameters for a process."""
        annotations = self.get_process_annotations(process_class)

        result = {
            'name': process_class.__name__,
            'description': annotations.get('summary', ''),
            'long_description': annotations.get('long', ''),
            'input': annotations.get('input', {}),
            'output': annotations.get('output', {}),
            'envs': annotations.get('envs', {})
        }

        return result
