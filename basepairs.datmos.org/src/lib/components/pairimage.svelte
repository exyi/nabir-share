<script lang="ts">
	import type { PairingInfo } from "$lib/pairing";
	import { getContext } from "svelte";
  import Modal, * as modal from 'svelte-simple-modal'
	import DetailModal from "./DetailModal.svelte";
  import type { Context } from 'svelte-simple-modal';


  export let url: string | undefined
  export let videoUrl: string | undefined
  export let pair: PairingInfo | undefined
  export let videoPreload: boolean = false
  export let videoShow: boolean = false
  export let allowHoverVideo: boolean = true
  export let parentSize: boolean = false
  export let linkText: string | undefined = undefined
  export let linkUrl: string | undefined = undefined
  export let onClick: () => any = () => { showModal(); return false }

  let pngFallback = false,
    webpFallback = false,
    videoLoaded = false
  const resolutions = [ [450, 800 ], [720, 1280], [1080, 1980], [1440, 2560 ] ],
    webpResolution = [ [450, 800 ], [1440, 2560 ] ]

  function generateSrcset(url: string, webpFallback: boolean) {
    const stripExt = url.replace(/\.\w+$/, '')
    const avifs = (webpFallback ? webpResolution : resolutions).map(r => `${stripExt}-${r[0]}.avif ${r[1]}w`).join(', ')
    return `${avifs}`
  }


  $: {
    url
    pngFallback = false
    webpFallback = false
    videoLoaded = false
  }

  function onerror(e: Event) {
    if (!webpFallback) {
      webpFallback = true
      return
    }
    if (!pngFallback) {
      pngFallback = true
      return
    }
    console.warn("Failed to load image", e)
  }

	function showModal(): any {
    const {open} = getContext<Context>('simple-modal');
    const imageUrl = pngFallback ? url : webpFallback ? url.replace(/\.\w+$/, '-1440.webp') : url.replace(/\.\w+$/, '-1440.avif')
    open(DetailModal, { pair, imageUrl, videoUrl }, {
      classContent: "smodal-content",
      styleWindow: {
        width: "80vw",
      }
    })
    return false
	}
</script>
<style>
  .pairimage {
    display: flex;
    flex-direction: column;
    border-radius: 15px;
  }
  .pairimage:hover {
    background-color: #dceff4;
    transition: background-color 0.2s;
    /* box-shadow: 0 0 5px 10px #dceff4; */
  }
  .pairimage.comparison-added {
    background-color: #e2f4dc;
    /* box-shadow: 0 0 5px 10px #e2f4dc; */
    box-shadow: inset 0 0 20px 10px white;
  }
  .pairimage.comparison-removed {
    background-color: #f4dce2;
    box-shadow: inset 0 0 10px 10px white;
  }
  .pairimage.comparison-added:hover {
    background-color: #c7e6c0;
    box-shadow: none;
  }
  .pairimage.comparison-removed:hover {
    background-color: #e6c0c7;
    box-shadow: none;
  }
  .header {
    display: flex;
    flex-direction: column;
    ;
    font-size: 0.75rem;
    margin-top: -1.5rem;
    text-align: center;
    height: 2rem;
  }
  .img-root {
    position: relative;
    cursor: pointer;
    aspect-ratio: 16 / 9;
  }
  @media (max-width: 800px) {
    .img-root.autosize {
      width: 49vw;
    }
  }
  @media (max-width: 1200px) and (min-width: 800px) {
    .img-root.autosize {
      width: 33vw;
    }
  }
  @media (max-width: 2000px) and (min-width: 1200px) {
    .img-root.autosize {
      width: 24vw;
    }
  }
  @media (min-width: 2000px) {
    .img-root.autosize {
      width: 19.5vw;
    }
  }
  .img-root .video {
    display: none;
  }
  .img-root.allow-video:hover .img, .img-root.video-show .img {
    display: none;
  }
  .img-root.allow-video:hover .video, .img-root.video-show .video {
    display: block;
  }
</style>

<div class="pairimage" class:comparison-added={pair.comparison === true} class:comparison-removed={pair.comparison === false}>
  <a on:click={() => onClick()} href={linkUrl ?? "javascript:;"}>
<div class="img-root" class:autosize={!parentSize} class:allow-video={allowHoverVideo && videoUrl != null} class:video-show={videoShow} on:mouseover={() => { videoLoaded = true }} on:focus={() => { videoLoaded = true }}>
  <div class="video">
    {#if videoLoaded || videoPreload || videoShow}
      <video src={videoUrl} autoplay loop muted preload="none"></video>
    {/if}
  </div>
  <div class="img">
    {#if url == null}
      <span>no image</span>
    {:else if pngFallback}
      <img src={url} alt="image wasn't generated" loading="lazy" />
    {:else}
      <img src={url} srcset={generateSrcset(url, webpFallback)} alt="xxx" loading="lazy" on:error={onerror} />
    {/if}
  </div>
</div>
</a>
<div class="header">
    <div style="flex: 1 1 0px;"></div>
    <div style="flex: 0 0 auto;">
      <strong>{pair?.id.nt1.pdbid}</strong>{pair?.id.nt1.model > 1 ? "-" + pair?.id.nt1.model : ""}
      {pair?.id.nt1.chain}-{pair?.id.nt1.resname??''}<strong>{pair?.id.nt1.resnum}{pair?.id.nt1.altloc??''}{pair?.id.nt1.inscode??''}</strong>
      · · ·
      {pair?.id.nt2.chain}-{pair?.id.nt2.resname??''}<strong>{pair?.id.nt2.resnum}{pair?.id.nt2.altloc??''}{pair?.id.nt2.inscode??''}</strong>
    </div>
    <div style="flex: 0 0 auto; overflow: visible">
      {#if linkText}
        <a href={linkUrl ?? "javascript:;"} style="text-decoration: underline;" on:click={() => onClick()}>{linkText}</a>
      {/if}
    </div>
</div>

</div>
