<script>
    import { Image, DataTable } from "$libs";
    import { Tabs, Tab, TabContent, Tile } from "$ccs";
</script>

{% for cdir in job.out.outdir | joinpaths: "*" | glob %}
<h1>{{cdir | basename}}</h1>

<h2>Chi-square test summary table</h2>

<DataTable
    src={{cdir | joinpaths: "chisq.txt" | quote}}
    data={ {{cdir | joinpaths: "chisq.txt" | datatable: sep="\t"}} }  />

<h2>Fisher's exact test summary table</h2>

<DataTable
    src={{cdir | joinpaths: "fisher.txt" | quote}}
    data={ {{cdir| joinpaths: "fisher.txt" | datatable: sep="\t"}} }  />


{%  for ddir in cdir | joinpaths: "*" | glob %}
<h2>{{ ddir | basename }}</h2>

{%    for edir in ddir | joinpaths: "*" | glob %}
{%      if ext(edir) == "txt" %}
{%          continue %}
{%      endif %}
<h3>{{ edir | basename }}</h3>

<h4> Contingency Table </h4>

<DataTable
    src={{edir | joinpaths: "contingency.txt" | quote}}
    data={ {{edir | joinpaths: "contingency.txt" | datatable: sep="\t"}} }  />

<h4> Bar Plot </h4>

<Image src={{edir | joinpaths: "counts.boxplot.png" | quote}} />

<h4> Statistic Tests </h4>

<Tabs>
    <Tab label="Chi-square Test" title="Chi-square Test" />
    <Tab label="Fisher's exact Test" title="Fisher's exact Test" />
    <div slot="content">
        <TabContent>
            <Tile><pre>{{ edir | joinpaths: "chisq.txt" | read | escape }}</pre></Tile>
        </TabContent>
        <TabContent>
            <Tile><pre>{{ edir | joinpaths: "fisher.txt" | read | escape }}</pre></Tile>
        </TabContent>
    </div>
</Tabs>

{%      endfor %}
{%  endfor %}
{% endfor %}
