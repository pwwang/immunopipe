# Immunopipe Agent Skills

## Overview

This directory contains **Agent Skills** for the immunopipe bioinformatics pipeline. These skills enable AI agents to generate immunopipe configuration files through natural language interactions.

**Total Skills**: 33 skills (1 main routing + 32 process-specific)  
**Format**: [Agent Skills specification](https://agentskills.io/specification) - directories with YAML frontmatter  
**Purpose**: AI-assisted TOML configuration generation

## Quick Start

### For AI Agents

1. Read `immunopipe-config/SKILL.md` - Main routing skill
2. Based on user requirements, route to appropriate process skill directories
3. Generate complete TOML configurations using skill documentation

### For Humans

Browse skill directories for comprehensive process documentation including:
- Complete configuration syntax
- All parameters with types, defaults, and descriptions
- Realistic usage examples (minimal, typical, advanced)
- Common usage patterns with complete TOML code
- Troubleshooting guides with solutions

## Skill Organization

```
skills/
├── immunopipe-config/             # Main routing skill (entry point)
│   └── SKILL.md                   # Pipeline architecture & process selection
├── sampleinfo/                    # Core workflow skills (8)
│   └── SKILL.md
├── seuratpreparing/
├── seuratclustering/
├── ...                            # + 29 more process skills
└── README.md                      # This file
```

**Note**: Each skill is a directory containing a `SKILL.md` file with YAML frontmatter following the [Agent Skills specification](https://agentskills.io/specification).

## Available Skills

### Core Workflow Skills (8)

Essential processes for basic scRNA-seq analysis:

| Skill Directory | Process | Purpose |
|-----------------|---------|---------|
| `sampleinfo/` | SampleInfo | Sample metadata entry point |
| `loadingrnafromseurat/` | LoadingRNAFromSeurat | Load pre-processed Seurat objects |
| `screploading/` | ScRepLoading | Load TCR/BCR repertoire data |
| `seuratpreparing/` | SeuratPreparing | QC, normalization, integration |
| `seuratclustering/` | SeuratClustering | Standard clustering workflow |
| `seuratclusteringofallcells/` | SeuratClusteringOfAllCells | Pre-selection clustering |
| `seuratsubclustering/` | SeuratSubClustering | Sub-clustering specific populations |
| `seuratclusterstats/` | SeuratClusterStats | Cluster visualization and QC |

### TCR/BCR Analysis Skills (7)

Processes for T cell and B cell receptor analysis:

| Skill Directory | Process | Purpose |
|-----------------|---------|---------|
| `torbcellselection/` | TOrBCellSelection | T or B cell selection from mixed populations |
| `screpcombiningexpression/` | ScRepCombiningExpression | Integrate TCR/BCR with RNA expression |
| `cdr3clustering/` | CDR3Clustering | CDR3 sequence similarity clustering |
| `tessa/` | TESSA | TCR epitope-specific sequence analysis |
| `cdr3aaphyschem/` | CDR3AAPhyschem | CDR3 physicochemical properties |
| `clonalstats/` | ClonalStats | Clonality and diversity metrics |

### Marker & Annotation Skills (7)

Processes for cell type identification and marker analysis:

| Skill Directory | Process | Purpose |
|-----------------|---------|---------|
| `clustermarkers/` | ClusterMarkers | Differential expression per cluster |
| `clustermarkersofallcells/` | ClusterMarkersOfAllCells | Markers before T/B selection |
| `celltypeannotation/` | CellTypeAnnotation | Automated cell type annotation |
| `seuratmap2ref/` | SeuratMap2Ref | Reference-based annotation |
| `markersfinder/` | MarkersFinder | Custom group comparisons |
| `topexpressinggenes/` | TopExpressingGenes | Highly expressed genes per cluster |
| `topexpressinggenesofallcells/` | TopExpressingGenesOfAllCells | Top genes before selection |

### Advanced Analysis Skills (5)

Specialized analysis processes:

| Skill Directory | Process | Purpose |
|-----------------|---------|---------|
| `modulescorecalculator/` | ModuleScoreCalculator | Gene signature scoring (exhaustion, proliferation, etc.) |
| `cellcellcommunication/` | CellCellCommunication | Ligand-receptor interaction inference |
| `cellcellcommunicationplots/` | CellCellCommunicationPlots | Communication network visualization |
| `scfgsea/` | ScFGSEA | Gene set enrichment analysis |
| `pseudobulkdeg/` | PseudoBulkDEG | Pseudo-bulk differential expression |

### Metabolic Analysis Skills (6)

Comprehensive metabolic pathway profiling:

| Skill Directory | Process | Purpose |
|-----------------|---------|---------|
| `scrnametaboliclandscape/` | ScrnaMetabolicLandscape | All-in-one metabolic analysis pipeline |
| `metabolicinput/` | MetabolicInput | Load metabolic pathway databases |
| `metabolicexpimputation/` | MetabolicExprImputation | Impute missing metabolic gene expression |
| `metabolicfeatures/` | MetabolicFeatures | Pathway enrichment analysis (GSEA-based) |
| `metabolicpathwayactivity/` | MetabolicPathwayActivity | Calculate pathway activity scores |
| `metabolicpathwayheterogeneity/` | MetabolicPathwayHeterogeneity | Analyze metabolic heterogeneity |

## Using Skills

### Option 1: MCP Server (Recommended)

```bash
# Start immunopipe MCP server
immunopipe mcp

# Configure in Claude Desktop or compatible AI agent
# Skills are automatically loaded and available
```

### Option 2: Direct File Access

AI agents with file system access can directly read skills:

```python
# Read main routing skill
main_skill = read("skills/immunopipe-config/SKILL.md")

# Read process-specific skills as needed
clustering_skill = read("skills/seuratclustering/SKILL.md")

# Generate configuration
config = generate_toml(user_request, main_skill, clustering_skill)
```

### Option 3: Manual Reference

Browse skills as comprehensive documentation:
- Each skill includes complete parameter reference
- Realistic configuration examples
- Common patterns and workflows
- Troubleshooting guidance

## Skill Structure

Each skill directory contains a `SKILL.md` file with:

**YAML Frontmatter** (required):
```yaml
---
name: skill-name  # Matches directory name
description: Brief description (max 1024 chars)
---
```

**Markdown Content**:
```markdown
# Process Name Configuration

## Purpose
What the process does and key capabilities

## When to Use
Conditions, requirements, and use cases

## Configuration Structure
Complete TOML syntax with all sections

## Configuration Examples
Minimal, typical, and advanced examples

## Common Patterns
Real-world usage scenarios with TOML code

## Dependencies
Upstream and downstream process relationships

## Validation Rules
Configuration constraints and requirements

## Troubleshooting
Common issues and solutions
```

## Example Workflows

### Workflow 1: Basic TCR-seq Analysis

**Goal**: Analyze PBMC TCR-seq data for clonal expansion

**Skills used**:
1. `immunopipe-config/` → Routes to TCR workflow
2. `sampleinfo/` → Configure samples
3. `seuratpreparing/` → QC and normalization
4. `screploading/` → Load TCR data
5. `seuratclustering/` → Cluster T cells
6. `screpcombiningexpression/` → Integrate TCR + RNA
7. `clonalstats/` → Clonality metrics

### Workflow 2: Exhausted T Cell Analysis

**Goal**: Identify and characterize exhausted CD8+ T cells

**Skills used**:
1. `seuratclustering/` → Cluster T cells
2. `modulescorecalculator/` → Score exhaustion signature
3. `clustermarkers/` → Find exhaustion-associated genes
4. `scfgsea/` → Pathway enrichment

### Workflow 3: Cell-Cell Communication

**Goal**: Analyze interactions between T cells and myeloid cells

**Skills used**:
1. `celltypeannotation/` → Annotate cell types
2. `cellcellcommunication/` → Infer interactions
3. `cellcellcommunicationplots/` → Visualize networks

## Skill Development

### Contributing

Skills are maintained in the immunopipe repository:
- **Location**: `skills/` directory (skill-per-directory structure)
- **Format**: [Agent Skills specification](https://agentskills.io/specification) - YAML frontmatter + Markdown
- **Updates**: Skills are updated when processes change
- **Pull requests**: Welcome for improvements and corrections

### Guidelines

When creating or updating skills:
1. **Directory structure**: Each skill must be in its own directory with a `SKILL.md` file
2. **YAML frontmatter**: Required `name` (matches directory) and `description` (max 1024 chars)
3. **Completeness**: Document all parameters with types and defaults
4. **Examples**: Include minimal, typical, and advanced configurations
5. **Clarity**: Use clear, concise language
6. **Accuracy**: Verify examples work with actual immunopipe
7. **Troubleshooting**: Include common issues from user experience

## Related Documentation

- **[Agent Skills Guide](../docs/agent-skills.md)**: Complete guide to using skills
- **[MCP Server](../docs/mcp-server.md)**: Model Context Protocol server setup
- **[Configurations](../docs/configurations.md)**: Manual TOML configuration reference
- **[Processes](../docs/processes/)**: Original process documentation

## External Resources

- [Agent Skills Specification](https://agentskills.io/specification) - Official format
- [Agent Skills Examples](https://github.com/anthropics/skills) - Anthropic examples
- [MCP Documentation](https://modelcontextprotocol.io) - Protocol specification

## Statistics

- **Total Skills**: 33 (1 main + 32 processes)
- **Total Lines**: ~12,700+ lines of comprehensive documentation
- **Coverage**: 100% of immunopipe processes (32/32)
- **Format**: [Agent Skills specification](https://agentskills.io/specification) compliant
- **Last Updated**: 2026-01-20

## Next Steps

1. **Browse skills**: Explore skill directories (e.g., `sampleinfo/`, `seuratclustering/`)
2. **Try MCP server**: Run `immunopipe mcp`
3. **Generate configs**: Use AI agents to create configurations
4. **Provide feedback**: Open issues for improvements