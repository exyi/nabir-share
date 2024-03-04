<script lang="ts">
	import { filterToSqlCondition, makeSqlQuery, type NucleotideFilterModel, type Range } from "$lib/dbModels";
  import type metadataModule from '$lib/metadata';
  import RangeSlider from 'svelte-range-slider-pips'

    export let filter: NucleotideFilterModel
    export let selectingFromTable: string | null = null
    export let metadata: typeof metadataModule[0] | null = null
    export let mode: "ranges" | "sql" = "ranges"

    let bonds = ["Bond 0", "Bond 1", "Bond 2"]
    $: bonds = metadata?.labels?.filter(l => l != null) ?? ["Bond 0", "Bond 1", "Bond 2"]
    function tryParseNum(x: string) {
      const n = Number.parseFloat(x)
      return isNaN(n) ? null : n
    }

    function ensureLength(array: Range[], index: number) {
      while (array.length <= index) {
        array.push({})
      }
    }
    // let ranges = {
    //     length: bonds.map(_ => ({ min: "", max: "" })),
    //     accAngle: bonds.map(_ => ({ min: "", max: "" })),
    //     donAngle: bonds.map(_ => ({ min: "", max: "" })),
    // }
    // function updateFilter() {
    //     function convertRange(x: any) {
    //         return { min: (isNaN(Number.parseFloat(x.min)) ? null : Number.parseFloat(x.min)),
    //                  max: (isNaN(Number.parseFloat(x.max)) ? null : Number.parseFloat(x.max)) }
    //     }
    //     filter = { ...filter,
    //         bond_length: bonds.map((_, i) => convertRange(ranges.length[i])),
    //         bond_acceptor_angle: bonds.map((_, i) => convertRange(ranges.accAngle[i])),
    //         bond_donor_angle: bonds.map((_, i) => convertRange(ranges.donAngle[i]))
    //     }
    // }

    // $:{
    //     ranges
    //     updateFilter()
    // }

    function modeChange(e: Event & { currentTarget:EventTarget & HTMLInputElement }) {
      mode = e.currentTarget.value as any

      const currentSqlQuery = makeSqlQuery(selectingFromTable.endsWith('_f') ? {... filter, filtered: false} : filter, selectingFromTable)
      if (mode=="sql" && !filter.sql) {
        filter = {...filter, sql: currentSqlQuery }
      }
      if (mode == "ranges" && filter.sql.trim() == currentSqlQuery.trim()) {
        filter = {...filter, sql: "" }
      }
    }

    function dnaRnaChange(e: Event & { currentTarget:EventTarget & HTMLInputElement }) {
      const v = {
        "rna": false,
        "dna": true,
        "both": undefined
      }[e.currentTarget.value as any]
      filter = {...filter, dna: v }
    }
</script>

<style>
  .panel-title {
    font-variant: small-caps;
    font-size: 1rem;
    font-weight: bold;
    text-align: center;
    height: 1.5rem;
  }
  .panel-field {
    height: 1.75rem;
    margin-bottom: 0.75rem;
  }
  .flex-columns {
    display: flex;
    width: 100%;
    justify-content: center;
    flex-wrap: wrap;
  }
  .flex-columns > * {
    flex: 0 1 auto;
  }

  .num-input {
    max-width: 100px;
  }

  .field-body {
    flex-grow: 1;
  }

  .mode-selection {
    display: flex;
    justify-content: center;
  }
  .sql-editor {
    font-family: 'Fira Code', 'Consolas', monospace;
  }
</style>

<div>
    <div class="control mode-selection">
        <label class="radio" title="Filter by constraining the H-bond parameters.">
          <input type="radio" checked={mode=="ranges"} value="ranges" name="editor_mode" on:change={modeChange}>
          Parameter ranges
        </label>
        <label class="radio" title="Filter by anything you want using SQL.">
          <input type="radio" checked={mode=="sql"} value="sql" name="editor_mode" on:change={modeChange}>
          SQL
        </label>
    </div>
    
    {#if mode=="ranges"}
    <div class="flex-columns" >
        <div class="column">
          <div class="panel-title"></div>
          {#each bonds as bond, i}
            <div class="panel-field">{bond}</div>
          {/each}
        </div>
        <div class="column">
            <h3 class="panel-title" title="Length in Å between the donor and acceptor heavy atoms">Length</h3>
            {#each bonds as bond, i}
              <!-- <div class="panel-field range-slider-wrapper">
                <RangeSlider min={0} max={6} step={0.1} pushy={true} suffix="Å" float={true}
                            values={[Number(ranges.length[i].min || 0), Number(ranges.length[i].max || 6)]}
                            on:change={e => { ranges.length[i].min = ""+e.detail.values[0]; ranges.length[i].max = ""+e.detail.values[1] }} />
              </div> -->
              <div class="panel-field field is-horizontal">
                <div class="field-body">
                  <div class="field">
                    <div class="field has-addons">
                      <p class="control">
                        <input class="input is-small num-input" type="number" step="0.1" min=0 max=6 placeholder="Min" value={filter.bond_length[i]?.min} on:change={ev => { ensureLength(filter.bond_length, i); filter.bond_length[i].min = tryParseNum(ev.currentTarget.value)} }>
                      </p>
                      <p class="control">
                        <input class="input is-small num-input" type="number" step="0.1" min=0 max=6 placeholder="Max" value={filter.bond_length[i]?.max} on:change={ev => { ensureLength(filter.bond_length, i); filter.bond_length[i].max = tryParseNum(ev.currentTarget.value)} }>
                      </p>
                    </div>
                  </div>
                </div>
              </div>
            {/each}
        </div>
        <div class="column">
            <h3 class="panel-title" title="Angle in degrees between the acceptor, donor and its covalently bound atom.">Donor Angle</h3>
            {#each bonds as bond, i}
            <div class="panel-field field is-horizontal">
                <!-- <div class="field-label">Bond {bond}</div> -->
                <div class="field-body">
                  <div class="field has-addons">
                    <p class="control">
                      <input class="input is-small num-input" type="number" step="5" min=0 max=360 placeholder="Min" value={filter.bond_donor_angle[i]?.min} on:change={ev => { ensureLength(filter.bond_donor_angle, i); filter.bond_donor_angle[i].min = tryParseNum(ev.currentTarget.value)} }>
                    </p>
                    <p class="control">
                      <input class="input is-small num-input" type="number" step="5" min=0 max=360 placeholder="Max" value={filter.bond_donor_angle[i]?.max} on:change={ev => { ensureLength(filter.bond_donor_angle, i); filter.bond_donor_angle[i].max = tryParseNum(ev.currentTarget.value)} }>
                    </p>
                  </div>
                </div>
              </div>
            {/each}
        </div>

        <div class="column">
            <h3 class="panel-title" title="Angle in degrees between the donor, acceptor and its covalently bound atom">Acceptor Angle</h3>
            {#each bonds as bond, i}
            <div class="panel-field field is-horizontal">
                <!-- <div class="field-label">Bond {bond}</div> -->
                <div class="field-body">
                  <div class="field has-addons">
                    <p class="control">
                      <input class="input is-small num-input" type="number" step="5" min=0 max=360 placeholder="Min" value={filter.bond_acceptor_angle[i]?.min} on:change={ev => { ensureLength(filter.bond_acceptor_angle, i); filter.bond_acceptor_angle[i].min = tryParseNum(ev.currentTarget.value)} }>
                    </p>
                    <p class="control">
                      <input class="input is-small num-input" type="number" step="5" min=0 max=360 placeholder="Max" value={filter.bond_acceptor_angle[i]?.max} on:change={ev => { ensureLength(filter.bond_acceptor_angle, i); filter.bond_acceptor_angle[i].max = tryParseNum(ev.currentTarget.value)} }>
                    </p>
                  </div>
                </div>
              </div>
            {/each}
        </div>

        <div class="column">
          <div class="control">
            <label class="radio" title="At least one of the nucleotides is RNA">
              <input type="radio" name="rna_dna_mode" value="rna" checked={filter.dna == false} on:change={dnaRnaChange}>
              RNA
            </label>
            <label class="radio" title="At least one of the nucleotides is DNA">
              <input type="radio" name="rna_dna_mode" value="dna" checked={filter.dna == true} on:change={dnaRnaChange}>
              DNA
            </label>
            <label class="radio">
              <input type="radio" name="rna_dna_mode" value="both" checked={filter.dna == null} on:change={dnaRnaChange}>
              Both
            </label>
          </div>
          <div class="field has-addons">
            {#if filter.resolution.min != null}
              <div class="control">
                <input class="input is-small num-input" style="max-width:4rem" type="number" step="0.1" min=0 max={filter.filtered ? 3.5 : 20} placeholder="Min" value={filter.resolution?.min ?? 0} on:change={ev => { filter.resolution ??= {}; filter.resolution.min = tryParseNum(ev.currentTarget.value)} }>
              </div>
            {/if}
            <label class="label" for="ntfilter-resolution">{#if filter.resolution.min != null}&nbsp;≤ {/if}Resolution ≤&nbsp;</label>
            <div class="control">
              <input class="input is-small num-input" style="max-width:4rem" type="number" step="0.1" min=0 max={filter.filtered ? 3.5 : 20} placeholder={filter.filtered ? '3.5' : ''} value={filter.resolution?.max ?? ''} on:change={ev => { filter.resolution ??= {}; filter.resolution.max = tryParseNum(ev.currentTarget.value)}}>
            </div>
            &nbsp;Å
          </div>

          <div class="control">
            <label class="checkbox" title="Filter out redundant nucleotides or nucleoties with bad something TODO">
              <input type="checkbox" checked={filter.filtered} on:change={e => filter = {...filter, filtered: e.currentTarget.checked }}>
              Representative set
            </label>

          {#if filter.filtered && filter.resolution.max && filter.resolution.max > 3.5}
            <p class="help is-danger">Representative set only<br> contains structures ≤3.5 Å</p>
          {/if}
          </div>
          <div class="control">
            <label class="checkbox" title="Include 'nearly pairs' as reported by fr3d basepair_detailed function">
              <input type="checkbox" checked={filter.includeNears} on:change={e => filter = {...filter, includeNears: e.currentTarget.checked }}>
              Include near-pairs
            </label>
          </div>
        </div>

        <div class="column">
          <div class="field">
            <label class="label" for="ntfilter-order-by">Order by</label>
            <div class="control">
              <div class="select is-small">
                <select bind:value={filter.orderBy} id="ntfilter-order-by">
                  <option value="">pdbid</option>
                  <option value="pdbid DESC, model DESC, chain1 DESC, nr1 DESC">pdbid descending</option>
                  <option value="resolution NULLS LAST, pdbid, model, chain1, nr1" title="Reported resolution of the source PDB structure">resolution</option>
                  <option value="mode_deviations" title="ASCENDING - best to worst - Number of standard deviations between the H-bond parameters and the modes (peaks) calculated from Kernel Density Estimate. Use to list &quot;nicest&quot; pairs and avoid secondary modes.">Deviation from KDE mode ↓</option>
                  <option value="mode_deviations" title="DESCENDING - worst to best - Number of standard deviations between the H-bond parameters and the modes (peaks) calculated from Kernel Density Estimate. Use to list &quot;nicest&quot; pairs and avoid secondary modes.">Deviation from KDE mode ↑</option>
                  <option value="log_likelihood DESC" title="↑ DESCENDING - best to worst - Multiplied likelihoods of all H-bond parameters in their Kernel Density Estimate distribution. Use to list &quot;nicest&quot; pairs without disqualifying secondary modes.">KDE likelihood ↑</option>
                  <option value="log_likelihood DESC" title="↓ ASCENDING - best to worst - Multiplied likelihoods of all H-bond parameters in their Kernel Density Estimate distribution. Use to list &quot;nicest&quot; pairs without disqualifying secondary modes.">KDE likelihood ↓</option>
                </select>
              </div>
            </div>
          </div>
        </div>

    </div>

    {:else if mode=="sql"}
      <div>
        <textarea class="textarea sql-editor" bind:value={filter.sql} style="width: 100%;"></textarea>
        <p class="help is-link">Use the SQL language to filter by anything. <code>selectedpair</code>, <code>selectedpair_f</code> and <code>selectedpair_n</code> contain the currently selected pair type, <code>_f</code> suffix are the filtered non-redundant set, <code>_n</code> suffix are the "nearly pairs". All other pair types are available in tables like <code>'tWW-A-A'</code> with the optional <code>_f</code> or <code>_n</code> suffix. The query runs in the browser, so run as many <code>DROP DATABASE</code>s as you please.</p>
        <p>
        </p>
      </div>
    {/if}
</div>
