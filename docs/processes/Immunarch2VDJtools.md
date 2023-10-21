# Immunarch2VDJtools

Convert immuarch format into VDJtools input formats.

This process converts the [`immunarch`](https://immunarch.com/) object to the
[`VDJtools`](https://vdjtools-doc.readthedocs.io/en/master/) input files,
in order to perform the VJ gene usage analysis by [`VJUsage`](./VJUsage.md) process.<br />

This process will generally generate a tab-delimited file for each sample,
with the following columns.<br />

- `count`: The number of reads for this clonotype
- `frequency`: The frequency of this clonotype
- `CDR3nt`: The nucleotide sequence of the CDR3 region
- `CDR3aa`: The amino acid sequence of the CDR3 region
- `V`: The V gene
- `D`: The D gene
- `J`: The J gene

See also: <https://vdjtools-doc.readthedocs.io/en/master/input.html#vdjtools-format>.<br />

This process has no environment variables.<br />

