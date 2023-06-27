{% from "utils/misc.liq" import report_jobs -%}

<script>
    import { Image, DataTable } from "$libs";
</script>


{%- macro report_job(job, h=1) -%}

{% for casedir in job.out.outdir | glob: "*" %}
<h{{h}}>{{casedir | stem}}</h{{h}}>

<Image src={{casedir | joinpaths: "overlapping_markers.png" | quote}} />

<DataTable
    src={{casedir | joinpaths: "overlapping_markers.txt" | quote }}
    data={ {{ casedir | joinpaths: "overlapping_markers.txt" | datatable: sep="\t", nrows=100 }} }
    />

{% endfor %}
{%- endmacro -%}


{%- macro head_job(job) -%}
<h1>Overlaping markers for multiple markers finder cases</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
