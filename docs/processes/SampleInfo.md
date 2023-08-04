# SampleInfo

This process is the entrance of the pipeline. It just pass by input file and list the sample information in the report.

To specify the input file:

```toml
[SampleInfo.in]
infile = [ "path/to/sample_info.txt" ]
```

Or with `pipen-board`, find the `SampleInfo` process and click the `Edit` button. Then you can specify the input file here:

![infile](images/SampleInfo-infile.png)

Theroetically, we can have multiple input files. However, it is not tested yet. If you have multiple input files to run, please run it with a different pipeline instance (configuration file).

For the content of the input file, please see details [here](../preparing-input.md#metadata).

Once the pipeline is finished, you can see the sample information in the report:

![report](images/SampleInfo-report.png)

Note that the required `RNAData` and `TCRData` columns are not shown in the report. They are used to specify the directories of the `scRNA-seq` and `scTCR-seq` data, respectively.
