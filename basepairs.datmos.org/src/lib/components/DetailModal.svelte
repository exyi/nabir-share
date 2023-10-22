<script lang="ts">
	import type { PairingInfo } from "$lib/pairing";

    export let imageUrl: string
    export let videoUrl: string
    export let pair: PairingInfo
</script>

<style>
    .imgpane {
        display: flex;
        align-items: center;
    }
    .imgpane > * {
        flex-grow: 1;
        width: 40vw;
    }


</style>
<div>
    <div class="imgpane">
        <img src={imageUrl} alt='x' />
        {#if videoUrl}
            <video src={videoUrl} autoplay muted loop />
        {/if}
    </div>
    <div>
        <pre>
            {JSON.stringify(pair, null, 2)}
        </pre>
    </div>
</div>
