<script>
    import { Image, DataTable } from "$libs";
</script>

<h1>Indicator gene expression vs clonotype percentage</h1>
{% for plotfile in job.out.outdir | glob: "*-vs-clonopct.png" %}
  <h2>{{ plotfile | basename | replace: "-vs-clonopct.png", "" }}</h2>
  <Image src={{plotfile | quote}} />
{% endfor %}

{%- if job.out.outdir | joinpaths: "kmeans.png" | as_path | attr: "is_file" | call -%}
  <h1>K-means clustering</h1>
  <p>When no `tcell_indicator` specified, a k-means clustering is done with `n=2`
    using the expressions of indicator genes and clonotype percentage for each cluster.
    Then the mean clonotype percentage is calculated for each kmeans cluster and the cluster with
    the highest mean clonotype percentage is selected, in which the Seurat clusters are
    assigned as T cells.
  </p>
  <Image src={{job.out.outdir | joinpaths: "kmeans.png" | quote}} />
{%- endif -%}

<h1>Data table</h1>

<DataTable data={ {{ job.out.outdir | joinpaths: "data.txt" | datatable: sep="\t" }} } />
