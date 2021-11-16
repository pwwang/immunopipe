If you have samples from the same patient, you can check the residency of the clonotypes (where the clones reside)

Clone residency is done by process `biopipen.namespaces.tcr.CloneResidency`. To check available configurations of it:

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
