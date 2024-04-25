<script lang="ts">
	import type { NumRange } from "$lib/dbModels";

  export let range: NumRange | null | undefined
  export let onChange: ((range: NumRange) => void) | undefined | null = null
  export let step = 1
  export let min: number | null = null
  export let max: number | null = null

  function tryParseNum(x: string) {
    const n = Number.parseFloat(x)
    return isNaN(n) ? null : n
  }

  function update(update: NumRange) {
    update = range ? { ...range, ...update } : update
    if (onChange) {
      onChange(update)
    } else {
      range = update
    }
  }

  function selectInput(ev: Event & { currentTarget: HTMLInputElement }) {
    ev.currentTarget.select()
  }
</script>

<style>
  .num-input {
    max-width: 100px;
  }
</style>

<div class="field is-horizontal">
  <div class="field-body">
    <div class="field">
      <div class="field has-addons">
      <p class="control">
          <input class="input is-small num-input" type="number" {step} {min} {max} placeholder="Min" value={range?.min} on:change={ev => update({ min: tryParseNum(ev.currentTarget.value) }) } on:click={selectInput}>
      </p>
      <p class="control">
          <input class="input is-small num-input" type="number" {step} {min} {max} placeholder="Max" value={range?.max} on:change={ev => update({ max: tryParseNum(ev.currentTarget.value) }) } on:click={selectInput}>
      </p>
      </div>
    </div>
  </div>
</div>
