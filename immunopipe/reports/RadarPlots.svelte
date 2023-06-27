<script>
    import { Image } from "$libs";
</script>

{% for figure in job.out.outdir | glob: "*.png" %}
    <div>
        <Image src={{figure | quote}} />
    </div>
{% endfor %}
