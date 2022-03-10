<script>
    import { Image, DataTable } from "@@";
    import { Tabs, Tab, TabContent, Tile } from "carbon-components-svelte";
</script>

<h1>Chi-square test summary table</h1>

<DataTable
    src={{job.out.outdir | joinpaths: "chisq.txt" | quote}}
    data={ {{job.out.outdir | joinpaths: "chisq.txt" | datatable: sep="\t"}} }  />

<h1>Fisher's exact test summary table</h1>

<DataTable
    src={{job.out.outdir | joinpaths: "fisher.txt" | quote}}
    data={ {{job.out.outdir | joinpaths: "fisher.txt" | datatable: sep="\t"}} }  />


{% for cdir in job.out.outdir | joinpaths: "*" | glob %}
<h1>{{cdir | basename}}</h1>

{%  for ddir in cdir | joinpaths: "*" | glob %}
<h2>{{ ddir | basename }}</h2>

<h3> Contingency Table </h3>

<DataTable
    src={{ddir | joinpaths: "contingency.txt" | quote}}
    data={ {{ddir | joinpaths: "contingency.txt" | datatable: sep="\t"}} }  />

<h3> Bar Plot </h3>

<Image src={{ddir | joinpaths: "counts.boxplot.png" | quote}} />

<h3> Statistic Tests </h3>

<Tabs>
    <Tab label="Chi-square Test" title="Chi-square Test" />
    <Tab label="Fisher's exact Test" title="Fisher's exact Test" />
    <div slot="content">
        <TabContent>
            <Tile><pre>{{ ddir | joinpaths: "chisq.txt" | read | escape }}</pre></Tile>
        </TabContent>
        <TabContent>
            <Tile><pre>{{ ddir | joinpaths: "fisher.txt" | read | escape }}</pre></Tile>
        </TabContent>
    </div>
</Tabs>

{%  endfor %}
{% endfor %}
