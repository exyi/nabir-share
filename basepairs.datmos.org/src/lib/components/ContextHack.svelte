<script lang="ts">
    import { getContext } from "svelte";

    export let name: string;
    export let value: any = undefined

    value = getContext(name)
</script>
