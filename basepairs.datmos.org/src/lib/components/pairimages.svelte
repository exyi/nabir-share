<script lang="ts">
    import Pairimage from "./pairimage.svelte";
    import type { PairId, PairingInfo } from "../pairing";
	import type { DetailModalViewModel } from "$lib/dbModels";

    export let pairs: PairingInfo[]
    export let rootImages: string
    export let rotImg: boolean = false
    export let imgAttachement: string
    export let videoAttachement: string
    export let videoOnHover: boolean = false
    export let onClick: (p: DetailModalViewModel) => void = () => {}

    function getUrl(pair: PairId, attachement:string, opt = { rotImg }) {
        if (!pair?.nt1.pdbid || !pair?.nt2.pdbid || pair.nt1.chain == null || pair.nt2.chain == null || pair.nt1.resnum == null || pair.nt2.resnum == null)
            return undefined
        const imgName = `${pair.nt1.chain}_${pair.nt1.resnum}${pair.nt1.inscode??''}${pair.nt1.altloc??''}-${pair.nt2.chain}_${pair.nt2.resnum}${pair.nt2.inscode??''}${pair.nt2.altloc??''}`
        return `${rootImages}/${pair.nt1.pdbid}/${imgName}${opt.rotImg ? "-rotX" : ""}${attachement}`
    }

</script>

<div class="imgcontainer">
    {#each pairs as p}
        <Pairimage pair={p} url={getUrl(p.id, imgAttachement)} videoUrl={getUrl(p.id, videoAttachement, {rotImg: false})} allowHoverVideo={videoOnHover}
            onClick={() => onClick({ pair: p, imgUrl: getUrl(p.id, imgAttachement, {rotImg: false}), rotImgUrl: getUrl(p.id, imgAttachement, {rotImg: true}), videoUrl: getUrl(p.id, videoAttachement, {rotImg: false}) }) } />
    {/each}
</div>

<style>
    .imgcontainer {
        display: flex;
        flex-wrap: wrap;
    }
</style>
