<script>
    import { Image } from "@@";
</script>

{% for job in jobs %}
<h1>{{ job.in.name }}</h1>

<Image src={{job.out.outfile | quote}} />
{% endfor %}
