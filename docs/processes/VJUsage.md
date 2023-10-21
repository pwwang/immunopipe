# VJUsage

Circos-style V-J usage plot displaying the frequency of various V-J junctions using vdjtools.

This process performs the VJ gene usage analysis using
[`VDJtools`](https://vdjtools-doc.readthedocs.io/en/master/).<br />
It wraps the [`PlotFancyVJUsage`](https://vdjtools-doc.readthedocs.io/en/master/basic.html#plotfancyvjusage) command in `VDJtools`.<br />
The output will be a V-J junction circos plot for a single sample.<br />
Arcs correspond to different V and J segments, scaled to their frequency in sample.<br />
Ribbons represent V-J pairings and their size is scaled to the pairing frequency
(weighted in present case).<br />

![VJUsage](https://vdjtools-doc.readthedocs.io/en/master/_images/basic-fancyvj.png){: width="80%" }

## Environment Variables

- `vdjtools`: *Default: `vdjtools`*. <br />
    The path to the `VDJtools` executable.<br />

