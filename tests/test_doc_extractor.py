"""Tests for documentation-based parameter extraction."""

from unittest.mock import Mock, patch
from immunopipe.mcp.doc_extractor import ProcessDocumentationExtractor

# Required for test ordering
pytest_order = 100


class TestProcessDocumentationExtractor:
    """Test the ProcessDocumentationExtractor class."""

    def test_init(self):
        """Test initialization."""
        extractor = ProcessDocumentationExtractor()
        assert extractor._process_annotations_cache == {}

    def test_extract_parameter_value_resolution(self):
        """Test extraction of resolution parameter."""
        extractor = ProcessDocumentationExtractor()

        # Test explicit resolution value
        result = extractor._extract_parameter_value(
            "resolution", {"type": "float"}, "clustering with resolution 0.5"
        )
        assert result == 0.5

        # Test high resolution suggestion
        result = extractor._extract_parameter_value(
            "resolution", {"type": "float"}, "high resolution clustering"
        )
        assert result == 1.2

        # Test low resolution suggestion
        result = extractor._extract_parameter_value(
            "resolution", {"type": "float"}, "low resolution clustering"
        )
        assert result == 0.3

    def test_extract_parameter_value_algorithm(self):
        """Test extraction of algorithm parameter."""
        extractor = ProcessDocumentationExtractor()

        # Test Leiden algorithm
        result = extractor._extract_parameter_value(
            "algorithm", {"type": "str"}, "use leiden algorithm for clustering"
        )
        assert result == "Leiden"

        # Test Louvain algorithm
        result = extractor._extract_parameter_value(
            "algorithm", {"type": "str"}, "use louvain algorithm for clustering"
        )
        assert result == "Louvain"

    def test_extract_parameter_value_cell_type(self):
        """Test extraction of cell_type parameter."""
        extractor = ProcessDocumentationExtractor()

        # Test T cell selection
        result = extractor._extract_parameter_value(
            "cell_type", {"type": "str"}, "select tcr data for analysis"
        )
        assert result == "T"

        # Test B cell selection
        result = extractor._extract_parameter_value(
            "cell_type", {"type": "str"}, "analyze bcr data"
        )
        assert result == "B"

    def test_convert_value(self):
        """Test value conversion."""
        extractor = ProcessDocumentationExtractor()

        # Test int conversion
        assert extractor._convert_value("42", "int") == 42
        assert extractor._convert_value("3.7", "int") == 3

        # Test float conversion
        assert extractor._convert_value("3.14", "float") == 3.14
        assert extractor._convert_value("42", "float") == 42.0

        # Test bool conversion
        assert extractor._convert_value("true", "bool") is True
        assert extractor._convert_value("false", "bool") is False
        assert extractor._convert_value("yes", "bool") is True
        assert extractor._convert_value("no", "bool") is False

        # Test string fallback
        assert extractor._convert_value("invalid_int", "int") == "invalid_int"

    @patch("immunopipe.mcp.doc_extractor.annotate")
    def test_get_process_annotations(self, mock_annotate):
        """Test getting process annotations."""
        extractor = ProcessDocumentationExtractor()

        # Mock process class
        mock_process = Mock()
        mock_process.__name__ = "TestProcess"

        # Mock annotate return value
        mock_annotate.return_value = {
            "envs": {
                "resolution": {"type": "float", "default": 0.8},
                "algorithm": {"type": "str", "default": "Leiden"},
            }
        }

        result = extractor.get_process_annotations(mock_process)

        assert "envs" in result
        assert "resolution" in result["envs"]
        assert "algorithm" in result["envs"]

        # Test caching
        result2 = extractor.get_process_annotations(mock_process)
        assert result == result2
        assert mock_annotate.call_count == 1  # Should be cached

    def test_extract_namespace_config(self):
        """Test extraction of namespace configuration."""
        extractor = ProcessDocumentationExtractor()

        # Use the actual pipen-annotate format
        ns_info = {
            "name": "FindClusters",
            "attrs": {"ns": True},
            "help": "Arguments for FindClusters resolution (type=float): "
            + "Clustering resolution. Default: 0.8 algorithm: Algorithm to use",
        }

        result = extractor._extract_namespace_config(
            "FindClusters",
            ns_info,
            "clustering with resolution 0.5 using leiden algorithm",
        )

        assert result is not None
        assert "resolution" in result
        assert result["resolution"] == 0.5
        assert "algorithm" in result
        assert result["algorithm"] == "Leiden"

    def test_extract_parameters_from_requirements_integration(self):
        """Test full integration of parameter extraction."""
        extractor = ProcessDocumentationExtractor()

        # Mock process class
        mock_process = Mock()
        mock_process.__name__ = "SeuratClustering"

        # Mock with actual pipen-annotate format (use Envs not envs)
        mock_annotations = {
            "Envs": {
                "FindClusters": {
                    "name": "FindClusters",
                    "attrs": {"ns": True},
                    "help": "Arguments for FindClusters resolution (type=float): "
                    + "Clustering resolution. Default: 0.8 algorithm: Algorithm to use",
                }
            }
        }

        extractor._process_annotations_cache[mock_process.__name__] = mock_annotations

        result = extractor.extract_parameters_from_requirements(
            "clustering with resolution 0.3 using louvain algorithm", mock_process
        )

        assert "FindClusters" in result
        assert result["FindClusters"]["resolution"] == 0.3
        assert result["FindClusters"]["algorithm"] == "Louvain"


def test_validation_function():
    """Validation function to test the doc extractor functionality."""
    from immunopipe.mcp.doc_extractor import ProcessDocumentationExtractor

    extractor = ProcessDocumentationExtractor()

    # Test basic functionality
    tests = [
        {
            "name": "resolution_extraction",
            "expected": 0.5,
            "actual": extractor._extract_parameter_value(
                "resolution", {"type": "float"}, "clustering with resolution 0.5"
            ),
        },
        {
            "name": "algorithm_extraction",
            "expected": "Leiden",
            "actual": extractor._extract_parameter_value(
                "algorithm", {"type": "str"}, "use leiden algorithm"
            ),
        },
        {
            "name": "cell_type_extraction",
            "expected": "T",
            "actual": extractor._extract_parameter_value(
                "cell_type", {"type": "str"}, "analyze tcr data"
            ),
        },
    ]

    failed_tests = []
    passed_tests = 0

    for test in tests:
        try:
            if test["actual"] == test["expected"]:
                passed_tests += 1
                print(f"✓ {test['name']}: PASSED")
            else:
                failed_tests.append(
                    f"{test['name']}: expected {test['expected']}, got {test['actual']}"
                )
                print(
                    f"✗ {test['name']}: FAILED - expected {test['expected']}, "
                    f"got {test['actual']}"
                )
        except Exception as e:
            failed_tests.append(f"{test['name']}: exception {str(e)}")
            print(f"✗ {test['name']}: FAILED - {str(e)}")

    total_tests = len(tests)
    print(f"\nTest Results: {passed_tests}/{total_tests} tests passed")

    if failed_tests:
        print("Failed tests:")
        for failure in failed_tests:
            print(f"  - {failure}")
        assert False

    print("All tests passed!")
    assert True


if __name__ == "__main__":
    test_validation_function()
