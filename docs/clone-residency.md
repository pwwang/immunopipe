If you have samples from the same patient, you can check the residency of the clonotypes (where the clones reside)

Clone residency is done by process [`biopipen.namespaces.tcr.CloneResidency`][1]. To check available configurations of it:

```shell
pipen run tcr CloneResidency
```

An example configuration could be:

```toml
# subject is specified at column "Patient" in metadata
subject = ["Patient"]
# where is the source of the clones
group = "Source"
# the order of the sources: X ~ Y
order = ["Normal", "Tumor"]
```

[1]: https://pwwang.github.io/biopipen/api/biopipen.ns.tcr/#biopipen.ns.tcr.CloneResidency
