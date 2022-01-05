{% from "utils/misc.liq" import report_jobs, table_of_images -%}

<script>
    import { Image, DataTable } from "@@";
    import { Tabs, Tab, TabContent } from "carbon-components-svelte";
</script>


{%- macro report_job(job, h=1) -%}
<h{{h}}>Markers intersection across comparisons</h{{h}}>
{% for casedir in job.out.outdir | joinpaths: "markers_overlapping", "*" | glob %}
<h{{h+1}}>{{ casedir | basename }}</h{{h+1}}>
<Tabs>
    <Tab label="Venn/Upset plot" title="Venn/Upset plot" />
    <Tab label="Data" title="Data" />
    <div slot="content">
        <TabContent>
            <Image src={{casedir | joinpaths: "overlapping_markers.png" | quote}} />
        </TabContent>
        <TabContent>
            <DataTable
                src={{casedir | joinpaths: "overlapping_markers.txt" | quote }}
                data={ {{ casedir | joinpaths: "overlapping_markers.txt" | datatable: sep="\t", nrows=100 }} }
                />
        </TabContent>
    </div>
</Tabs>
{% endfor %}

<h{{h}}>Meta markers of &gt; 2 groups by ANOVA</h{{h}}>
{% for casedir in job.out.outdir | joinpaths: "meta_markers", "*" | glob %}
<h{{h+1}}>{{casedir | basename}}</h{{h+1}}>
<DataTable
    src={{casedir | joinpaths: "meta_makers.txt" | quote }}
    data={ {{ casedir | joinpaths: "meta_makers.txt" | datatable: sep="\t", nrows=100 }} }
    />

{% set images = [] %}
{% set names = [] %}

{% for imgpath in casedir | joinpaths: "*.png" | glob %}
{% set _ = images.append(imgpath) %}
{% set _ = names.append(stem(imgpath)) %}
{% endfor %}

{{ table_of_images(images, names) }}
{% endfor %}
{%- endmacro -%}


{%- macro head_job(job) -%}
{% assign config = job.in.configfile | read | toml_loads %}
<h1>{{config.name | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
