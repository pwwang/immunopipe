# Immunarch2VDJtools

This process converts the [`immunarch`][1] object to the [`VDJtools`][2] input files, in order to perform the VJ gene usage analysis by [`VJUsage`][3] process.

This process will generally generate a tab-delimited file for each sample, with the following columns:

- `count`: The number of reads for this clonotype
- `frequency`: The frequency of this clonotype
- `CDR3nt`: The nucleotide sequence of the CDR3 region
- `CDR3aa`: The amino acid sequence of the CDR3 region
- `V`: The V gene
- `D`: The D gene
- `J`: The J gene

See also: <https://vdjtools-doc.readthedocs.io/en/master/input.html#vdjtools-format>.

This process has no environment variables.


[1]: https://immunarch.com/
[2]: https://vdjtools-doc.readthedocs.io/en/master/
[3]: ./VJUsage.md
