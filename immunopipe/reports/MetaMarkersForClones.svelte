{% from "utils/misc.liq" import report_jobs, table_of_images -%}

<script>
    import { Image, DataTable } from "@@";
</script>


{%- macro report_job(job, h=1) -%}

{% for casedir in job.out.outdir | glob: "*" %}
<h{{h}}>{{casedir | joinpaths: "case.txt" | read | escape}}</h{{h}}>
<DataTable
    src={{casedir | joinpaths: "meta_makers.txt" | quote }}
    data={ {{ casedir | joinpaths: "meta_makers.txt" | datatable: sep="\t", nrows=100 }} }
    />

{{ table_of_images(glob(joinpaths(casedir, "*.png"))) }}
{% endfor %}
{%- endmacro -%}


{%- macro head_job(job) -%}
<h1>Meta-markers for multiple groups</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
