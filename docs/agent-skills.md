# Agent Skills

## Overview

**Agent Skills** are AI-powered configuration assistants that help you create immunopipe configuration files through natural language. Instead of manually writing TOML configuration files, you can describe what you want to do in plain English, and AI agents use these skills to generate the correct configuration.

/// Note | What are Agent Skills?
Agent Skills are a lightweight, open format for extending AI agent capabilities with specialized knowledge and workflows. Learn more at [agentskills.io](https://agentskills.io).
///

## Why Use Agent Skills?

**Problem**: Immunopipe has 30+ processes with hundreds of configuration options. Understanding which processes to use, how to configure them, and how they interact requires deep pipeline knowledge.

**Solution**: Agent Skills provide structured knowledge that AI agents can use to:
- **Recommend processes** based on your analysis goals
- **Generate complete configurations** with sensible defaults
- **Explain parameters** and their effects
- **Validate configurations** before running

## Available Skills

Immunopipe provides two types of skills:

### 1. Main Pipeline Skill

**Location**: `skills/immunopipe-config/SKILL.md`

This is the entry point skill that helps AI agents:
- Assess your data type (scRNA-seq only vs scRNA+scTCR/BCR-seq)
- Determine analysis goals (QC, clustering, TCR analysis, metabolic profiling, etc.)
- Route to appropriate process-specific skills
- Generate pipeline-level configuration

**Example use**: "I have PBMC data with TCR-seq. I want to identify exhausted CD8+ T cells and analyze their clonality."

### 2. Process-Specific Skills

**Location**: `skills/{process-name}/SKILL.md` (32 individual process skills)

Each process has a dedicated skill documenting:
- **Purpose**: What the process does
- **When to use**: Conditions and requirements
- **Configuration structure**: Complete TOML syntax
- **Parameters**: All environment variables with types, defaults, and descriptions
- **Examples**: Minimal, typical, and advanced configurations
- **Common patterns**: Real-world usage scenarios
- **Troubleshooting**: Common issues and solutions

**Categories**:

#### Core Workflow (8 processes)
- `SampleInfo` - Sample metadata entry point
- `LoadingRNAFromSeurat` - Load pre-processed Seurat objects
- `ScRepLoading` - Load TCR/BCR repertoire data
- `SeuratPreparing` - QC, normalization, integration
- `SeuratClustering` - Standard clustering workflow
- `SeuratClusteringOfAllCells` - Pre-selection clustering
- `SeuratSubClustering` - Sub-clustering specific populations
- `SeuratClusterStats` - Cluster visualization and QC

#### TCR/BCR Analysis (7 processes)
- `TOrBCellSelection` - T or B cell selection
- `ScRepCombiningExpression` - Integrate TCR/BCR with RNA
- `CDR3Clustering` - CDR3 sequence clustering
- `TESSA` - TCR epitope-specific analysis
- `CDR3AAPhyschem` - Physicochemical properties
- `ClonalStats` - Clonality and diversity metrics

#### Marker & Annotation (7 processes)
- `ClusterMarkers` - Differential expression per cluster
- `ClusterMarkersOfAllCells` - Markers before selection
- `CellTypeAnnotation` - Automated cell type annotation
- `SeuratMap2Ref` - Reference-based annotation
- `MarkersFinder` - Custom group comparisons
- `TopExpressingGenes` - Highly expressed genes
- `TopExpressingGenesOfAllCells` - Top genes before selection

#### Advanced Analysis (10 processes)
- `ModuleScoreCalculator` - Gene signature scoring
- `CellCellCommunication` - Ligand-receptor interactions
- `CellCellCommunicationPlots` - Communication visualization
- `ScFGSEA` - Gene set enrichment analysis
- `PseudoBulkDEG` - Pseudo-bulk differential expression
- Metabolic processes (6): `ScrnaMetabolicLandscape`, `MetabolicInput`, `MetabolicExprImputation`, `MetabolicFeatures`, `MetabolicPathwayActivity`, `MetabolicPathwayHeterogeneity`

## How to Use Agent Skills

### Option 1: MCP Server (Recommended for Claude Desktop)

Immunopipe includes a built-in Model Context Protocol (MCP) server that exposes skills to AI agents like Claude:

```bash
# Start MCP server
immunopipe mcp

# Or configure in Claude Desktop settings
# See docs/mcp-server.md for full setup
```

The MCP server provides tools:
- `configure_immunopipe` - Generate configurations via natural language
- `generate_full_config` - Create complete pipeline configurations
- `suggest_processes` - Recommend processes for your analysis
- `list_options` - Browse available configuration options

**Example interaction**:
```
User: I have 5 PBMC samples with TCR-seq data. I want to:
1. QC and normalize the data
2. Cluster T cells
3. Identify exhausted T cells
4. Analyze TCR clonality by cluster

Agent: I'll generate a configuration for your TCR-seq analysis workflow.
      Let me use the immunopipe configuration skill...

[Agent generates complete TOML configuration with:
 - SampleInfo for your 5 samples
 - SeuratPreparing with appropriate QC thresholds
 - SeuratClustering for T cell clustering
 - ModuleScoreCalculator for exhaustion scoring
 - ScRepCombiningExpression to integrate TCR data
 - ClonalStats for clonality analysis]

Here's your configuration. Would you like me to explain any section?
```

### Option 2: Direct Skill Access (For Other AI Agents)

AI agents with file system access can directly read skill files:

```python
# Example: Agent reading skills
main_skill = read("skills/immunopipe-config/SKILL.md")
process_skill = read("skills/seuratclustering/SKILL.md")

# Agent uses skill content to generate configuration
config = generate_toml_from_user_request(user_query, main_skill, process_skill)
```

### Option 3: Manual Reference

Even without AI agents, skills serve as comprehensive documentation:
- Browse skill directories (e.g., `skills/seuratclustering/`, `skills/clustermarkers/`)
- Each `SKILL.md` file includes complete configuration examples
- Copy-paste examples and adapt to your needs

## Skill Format

Skills follow the [Agent Skills specification](https://agentskills.io/specification):

**Directory Structure**:
```
skills/
├── immunopipe-config/
│   └── SKILL.md                # Main routing skill
├── sampleinfo/
│   └── SKILL.md                # Process-specific skill
├── seuratclustering/
│   └── SKILL.md
└── ...                         # + 30 more process skills
```

**File Format** (`SKILL.md`):
```markdown
---
name: process-name              # Matches directory name
description: When to use this process and what it does (max 1024 chars)
---

# Process Name

## Purpose
What this process does...

## Configuration Examples
[Complete TOML examples]

## Troubleshooting
[Common issues and solutions]
```

## Example Workflows

### Workflow 1: Basic TCR-seq Analysis

**User goal**: "Analyze PBMC TCR-seq data to identify clonally expanded T cells"

**Skills used**:
1. `immunopipe-config/` - Routes to TCR workflow
2. `sampleinfo/` - Configure sample metadata
3. `seuratpreparing/` - QC and normalization
4. `screploading/` - Load TCR data
5. `seuratclustering/` - Cluster T cells
6. `screpcombiningexpression/` - Integrate TCR + RNA
7. `clonalstats/` - Calculate clonality metrics

**Output**: Complete TOML configuration with all processes correctly wired

### Workflow 2: Cell-Cell Communication

**User goal**: "Analyze ligand-receptor interactions between T cells and myeloid cells"

**Skills used**:
1. `immunopipe-config/` - Routes to non-TCR workflow
2. `seuratclustering/` - Cluster all cells
3. `celltypeannotation/` - Annotate T and myeloid cells
4. `cellcellcommunication/` - Infer interactions
5. `cellcellcommunicationplots/` - Visualize networks

**Output**: Configuration focused on communication analysis

### Workflow 3: Metabolic Profiling

**User goal**: "Compare metabolic pathway activity between exhausted and effector T cells"

**Skills used**:
1. `modulescorecalculator/` - Score exhaustion/effector signatures
2. `scrnametaboliclandscape/` - Full metabolic analysis pipeline
3. `seuratclusterstats/` - Visualize metabolic scores

**Output**: Configuration with metabolic analysis modules

## Benefits for Different Users

### For Beginners
- **Lower barrier**: Describe goals in natural language instead of learning TOML syntax
- **Guided selection**: AI recommends appropriate processes
- **Error prevention**: Generated configs are validated
- **Learning resource**: Skills document all options with examples

### For Experienced Users
- **Faster prototyping**: Generate base configs quickly, then customize
- **Parameter discovery**: Find advanced options you didn't know existed
- **Best practices**: Learn common patterns from skill examples
- **Troubleshooting**: Quick reference for common issues

### For Workflow Developers
- **Reproducibility**: AI-generated configs are complete and documented
- **Standardization**: Consistent configuration patterns across projects
- **Collaboration**: Share natural language requirements, generate consistent configs

## Skill Development

Immunopipe skills are actively maintained alongside the pipeline:

- **Location**: [`skills/`](https://github.com/pwwang/immunopipe/tree/dev/skills) directory in repository
- **Format**: [Agent Skills specification](https://agentskills.io/specification) - directories with `SKILL.md` files
- **Structure**: One directory per skill with YAML frontmatter
- **Updates**: Skills are updated when processes change
- **Contributions**: Submit improvements via pull requests

Each skill includes:
- YAML frontmatter with `name` and `description`
- Complete process documentation
- All configuration parameters with types and defaults
- Realistic usage examples (minimal, typical, advanced)
- Common troubleshooting scenarios with solutions

## Related Documentation

- **[MCP Server](mcp-server.md)**: Set up the Model Context Protocol server for Claude Desktop
- **[Configurations](configurations.md)**: Manual TOML configuration reference
- **[Processes](processes/SampleInfo.md)**: Individual process documentation
- **[Getting Started](getting-started.md)**: Quick start guide for immunopipe

## External Resources

- [Agent Skills Specification](https://agentskills.io/specification) - Official format documentation
- [Agent Skills Examples](https://github.com/anthropics/skills) - Example skills from Anthropic
- [MCP Documentation](https://modelcontextprotocol.io) - Model Context Protocol specification