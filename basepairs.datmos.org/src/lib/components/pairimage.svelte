<script lang="ts">
	import type { PairingInfo } from "$lib/pairing";

  export let url: string | undefined
  export let videoUrl: string | undefined
  export let pair: PairingInfo | undefined

  let imgResolution
</script>
<style>
  .pairimage {
    display: flex;
    flex-direction: column;
    flex-grow: 1;
  }
  .img-root {
    position: relative;
    aspect-ratio: 16 / 9;
  }
  @media (max-width: 800px) {
    .img-root {
      width: 50vw;
    }
  }
  @media (max-width: 1200px) and (min-width: 800px) {
    .img-root {
      width: 33vw;
    }
  }
  @media (max-width: 2000px) and (min-width: 1200px) {
    .img-root {
      width: 25vw;
    }
  }
  @media (min-width: 2000px) {
    .img-root {
      width: 20vw;
    }
  }
  .img-root .video {
    display: none;
  }
  .img-root:hover .img {
    /* opacity: 0; */
    display: none;
  }
  .img-root:hover .video {
    display: block;
    /* width: 0;
    height: 0;
    position: absolute;
    top: 0;
    left: 0;
    overflow: visible; */
  }
</style>

<div class="pairimage">
<div class="header">
    {pair?.id.nt1.pdbid}{pair?.id.nt1.model > 1 ? "-" + pair?.id.nt1.model : ""}
    {pair?.id.nt1.chain}-{pair?.id.nt1.resname??''}{pair?.id.nt1.resnum}{pair?.id.nt1.altloc??''}{pair?.id.nt1.inscode??''}
    · · ·
    {pair?.id.nt2.chain}-{pair?.id.nt2.resname??''}{pair?.id.nt2.resnum}{pair?.id.nt2.altloc??''}{pair?.id.nt2.inscode??''}
</div>
<div class="img-root">
  <div class="video">
    <video src={videoUrl} autoplay loop muted></video>
  </div>
  <div class="img">
    <img src={url} alt="x" loading="lazy" />
  </div>
</div>

</div>
