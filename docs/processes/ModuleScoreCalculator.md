# ModuleScoreCalculator

/// Tip | `ModuleScoreCalculator` is added in `0.7.0`
///

This process calculates the module scores for each cell

The module scores are calculated by [`Seurat::AddModuleScore()`](https://satijalab.org/seurat/reference/addmodulescore) or [`Seurat::CellCycleScoring()`](https://satijalab.org/seurat/reference/cellcyclescoring) for cell cycle scores.

The module scores are calculated as the average expression levels of each program on single cell level, subtracted by the aggregated expression of control feature sets. All analyzed features are binned based on averaged expression, and the control features are randomly selected from each bin.

## Environment Variables

- `defaults` (`ns`): The default parameters for modules.
    - `features`: The features to calculate the scores. Multiple features should be separated by comma. You can also specify cc.genes or cc.genes.updated.2019 to use the cell cycle genes to calculate cell cycle scores. If so, three columns will be added to the metadata, including S.Score, G2M.Score and Phase. Only one type of cell cycle scores can be calculated at a time.
    - `nbin` (`type=int`): Number of bins of aggregate expression levels for all analyzed features.
    - `ctrl` (`type=int`): Number of control features selected from the same bin per analyzed feature.
    - `k` (`flag`): Use feature clusters returned from DoKMeans.
    - `assay`: The assay to use.
    - `seed` (`type=int`): Set a random seed.
    - `search` (`flag`): Search for symbol synonyms for features in features that don't match features in object?
    - `keep` (`flag`): Keep the scores for each feature? Only works for non-cell cycle scores.
    - `agg` (`choice`): The aggregation function to use. Only works for non-cell cycle scores.
        - `mean`: The mean of the expression levels
        - `median`: The median of the expression levels
        - `sum`: The sum of the expression levels
        - `max`: The max of the expression levels
        - `min`: The min of the expression levels
        - `var`: The variance of the expression levels
        - `sd`: The standard deviation of the expression levels
- `modules` (`type=json`): The modules to calculate the scores. Keys are the names of the expression programs and values are the dicts inherited from env.defaults. Here are some examples -
    ```json
    {
        "CellCycle": {"features": "cc.genes.updated.2019"},
        "Exhaustion": {"features": "HAVCR2,ENTPD1,LAYN,LAG3"},
        "Activation": {"features": "IFNG"},
        "Proliferation": {"features": "STMN1,TUBB"}
    }
    ```
