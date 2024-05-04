<script lang="ts">
	import type { NumRange } from "$lib/dbModels";

  export let clear1: NumRange | null | undefined = null
  export let clear2: NumRange | null | undefined = null
  export let clear3: NumRange | null | undefined = null
  export let clearList: NumRange[] | null | undefined | any = null
  export let hide: boolean = true

  let any = false
  $: any = nn(clear1) || nn(clear2) || nn(clear3) || clearList?.some(nn)

  function nn(v: NumRange | null | undefined): v is NumRange {
    return v && (v.min != null || v.max != null)
  }
</script>

<style>
  button:disabled {
    filter: grayscale(1);
  }
  button {
    background: none;
    border: none;
    cursor: pointer;
    padding: 0;
    margin-left: 0.5rem;
  }
</style>

<span style="float: right">
  {#if !hide || any}
  <button
    disabled={!any}
    on:click={() => {
      clear1 = {}
      clear2 = {}
      clear3 = {}
      if (clearList)
        clearList = []
    }}
    >
    ❌
  </button>
  {/if}
</span>
