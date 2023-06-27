<script>
    import { DataTable } from "$libs";
</script>

<h1>Sample Information</h1>

<DataTable
    data={ {{ job.out.outfile | datatable: sep="\t", excluded=["TCRDir", "RNADir"] }} }
    pageSize={50}
    />


