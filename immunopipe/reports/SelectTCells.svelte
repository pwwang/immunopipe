<script>
    import { Image, DataTable } from "@@";
</script>

<h1>Division rational</h1>

<Image src={{job.out.outdir | joinpaths: "subdivision-rationale.png" | quote}} />

<DataTable data={ {{ job.out.outdir | joinpaths: "clusters.txt" | datatable: sep="\t" }} } />
