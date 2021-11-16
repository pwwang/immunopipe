<script>
    import { DataTable } from "@@";
</script>

<h1>Sample Information</h1>

<DataTable
    data={ {{ job.out.outfile | datatable: sep="\t", ncols=[0,1,2,5,6,7,8,9,10,11,12] }} }
    pageSize={50}
    />


