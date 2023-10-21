# ModuleScoreCalculator

Calculate the module scores for each cell

/// Tip | `ModuleScoreCalculator` is added in `0.7.0`
///

The module scores are calculated by
[`Seurat::AddModuleScore()`](https://satijalab.org/seurat/reference/addmodulescore)
or [`Seurat::CellCycleScoring()`](https://satijalab.org/seurat/reference/cellcyclescoring)
for cell cycle scores.<br />

The module scores are calculated as the average expression levels of each
program on single cell level, subtracted by the aggregated expression of
control feature sets. All analyzed features are binned based on averaged
expression, and the control features are randomly selected from each bin.<br />

## Environment Variables

- `defaults` *(`ns`)*:
    The default parameters for `modules`.<br />
    - `features`:
        The features to calculate the scores. Multiple features
        should be separated by comma.<br />
        You can also specify `cc.genes` or `cc.genes.updated.2019` to
        use the cell cycle genes to calculate cell cycle scores.<br />
        If so, three columns will be added to the metadata, including
        `S.Score`, `G2M.Score` and `Phase`.<br />
        Only one type of cell cycle scores can be calculated at a time.<br />
    - `nbin` *(`type=int`)*: *Default: `24`*. <br />
        Number of bins of aggregate expression levels
        for all analyzed features.<br />
    - `ctrl` *(`type=int`)*: *Default: `100`*. <br />
        Number of control features selected from
        the same bin per analyzed feature.<br />
    - `k` *(`flag`)*: *Default: `False`*. <br />
        Use feature clusters returned from `DoKMeans`.<br />
    - `assay`: *Default: `RNA`*. <br />
        The assay to use.<br />
    - `seed` *(`type=int`)*: *Default: `8525`*. <br />
        Set a random seed.<br />
    - `search` *(`flag`)*: *Default: `False`*. <br />
        Search for symbol synonyms for features in
        features that don't match features in object?<br />
    - `keep` *(`flag`)*: *Default: `False`*. <br />
        Keep the scores for each feature?<br />
        Only works for non-cell cycle scores.<br />
    - `agg` *(`choice`)*: *Default: `mean`*. <br />
        The aggregation function to use.<br />
        Only works for non-cell cycle scores.<br />
        - `mean`:
            The mean of the expression levels
        - `median`:
            The median of the expression levels
        - `sum`:
            The sum of the expression levels
        - `max`:
            The max of the expression levels
        - `min`:
            The min of the expression levels
        - `var`:
            The variance of the expression levels
        - `sd`:
            The standard deviation of the expression levels
- `modules` *(`type=json`)*: *Default: `{}`*. <br />
    The modules to calculate the scores.<br />
    Keys are the names of the expression programs and values are the
    dicts inherited from `env.defaults`.<br />
    Here are some examples -

    ```python
    {
        "CellCycle": {"features": "cc.genes.updated.2019"},
        "Exhaustion": {"features": "HAVCR2,ENTPD1,LAYN,LAG3"},
        "Activation": {"features": "IFNG"},
        "Proliferation": {"features": "STMN1,TUBB"}
    }
    ```


