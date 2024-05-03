<script lang="ts">
	import type { PairId } from '$lib/pairing';
    import { onMount } from 'svelte';

    let container: HTMLDivElement;
    export let pairId: PairId;

    onMount(() => {
        import('./molStarBackend').then(({ createMolStar }) => {
            createMolStar(container, pairId);
        });
    });

</script>


<div bind:this={container}></div>
