# VJUsage

This process performs the VJ gene usage analysis using [`VDJtools`][1]. It wraps the [`PlotFancyVJUsage`][2] command in `VDJtools`. The output will be a V-J junction circos plot for a single sample. Arcs correspond to different V and J segments, scaled to their frequency in sample. Ribbons represent V-J pairings and their size is scaled to the pairing frequency (weighted in present case).

![VJUsage](https://vdjtools-doc.readthedocs.io/en/master/_images/basic-fancyvj.png){: width="80%" }

## Environment variables

- `vdjtools`: The path to the `VDJtools` executable. Default: `vdjtools`.
- `vdjtools_patch`: The patch file for `VDJtools`. It's delivered with the pipeline ([`biopipen`][3] package).
    - You don't need to provide this file, unless you want to use a different patch file by yourself.
    - See the issue with `VDJtools` [here](https://github.com/mikessh/vdjtools/issues/139).


[1]: https://vdjtools-doc.readthedocs.io/en/master/
[2]: https://vdjtools-doc.readthedocs.io/en/master/basic.html#plotfancyvjusage
[3]: https://github.com/pwwang/biopipen
