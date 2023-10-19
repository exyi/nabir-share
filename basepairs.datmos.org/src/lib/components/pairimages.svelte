<script lang="ts">
    import Pairimage from "./pairimage.svelte";
    import type { PairId, PairingInfo, PairingType } from "../pairing";

    export let pairs: PairingInfo[]
    export let rootImages: string
    export let imgAttachement: string
    export let videoAttachement: string

    function getUrl(pair: PairId, attachement:string) {
        const imgName = `${pair.nt1.chain}_${pair.nt1.resnum}${pair.nt1.inscode??''}${pair.nt1.altloc??''}-${pair.nt2.chain}_${pair.nt2.resnum}${pair.nt2.inscode??''}${pair.nt2.altloc??''}`
        return `${rootImages}/${pair.nt1.pdbid}/${imgName}${attachement}`
    }

</script>

<div class="imgcontainer">
    {#each pairs as p}
        <Pairimage pair={p} url={getUrl(p.id, imgAttachement)} videoUrl={getUrl(p.id, videoAttachement)} />
    {/each}
</div>

<style>
    .imgcontainer {
        display: flex;
        flex-wrap: wrap;
    }
</style>
