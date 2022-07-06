Basic TCR analysis uses [immunarch][1] to analyze TCR data.

The analyses include:

- Basic analysis and clonality
- Repertoire overlap and public clonotypes
- Gene usage analysis
- Diversity estimation
- Tracking of clonotypes
- Kmer and sequence motif analysis and visualisation

The analyses are wrapped with a `pipen` process. See also: [biopipen.ns.tcr.Immunarch][2]

You could also run
```shell
pipen run tcr Immunarch
```
to check configurable items.

[1]: https://immunarch.com/
[2]: https://pwwang.github.io/biopipen/api/biopipen.ns.tcr/#biopipen.ns.tcr.Immunarch
