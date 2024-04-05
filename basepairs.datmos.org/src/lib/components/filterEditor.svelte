<script lang="ts">
	import { defaultFilter, filterToSqlCondition, getDataSourceTable, makeSqlQuery, type ComparisonMode, type NucleotideFilterModel, type NumRange, orderByOptions } from "$lib/dbModels";
  import type metadataModule from '$lib/metadata';
  import RangeSlider from 'svelte-range-slider-pips'
  import * as filterfilter from '$lib/predefinedFilterLoader'
	import RangeEditor from "./RangeEditor.svelte";
	import _ from "lodash";

    export let filter: NucleotideFilterModel
    export let filterBaseline: NucleotideFilterModel | undefined
    export let allowFilterBaseline: boolean = true
    export let comparisonMode: ComparisonMode | undefined = undefined
    export let selectingFromTable: string | null = null
    export let metadata: typeof metadataModule[0] | null = null
    export let mode: "ranges" | "sql" | "basic" = "basic"

    let bonds = ["Bond 0", "Bond 1", "Bond 2"]
    $: bonds = metadata?.labels?.filter(l => l != null) ?? ["Bond 0", "Bond 1", "Bond 2"]
    function tryParseNum(x: string) {
      const n = Number.parseFloat(x)
      return isNaN(n) ? null : n
    }

    function setBaseline(f: NucleotideFilterModel | undefined, mode: string) {
      if (mode == "sql" || f == null) {
        filterBaseline = f
      } else {
        filterBaseline = { ...f, sql: undefined }
      }
    }

    function formatFilter(f: NucleotideFilterModel, compareWith: NucleotideFilterModel | undefined) {
      if (f == null) return ""
      const clauses = filterToSqlCondition(f).filter(x => !["jirka_approves"].includes(x))
      const datasetName = f.datasource == "fr3d-f" ? "FR3D, Representative Set" :
                          f.datasource == "fr3d" ? "3D, entire PDB" :
                          f.datasource == "fr3d-nf" ? "FR3D with nears, RS" :
                          f.datasource == "fr3d-n" ? "FR3D with nears, PDB" :
                          f.datasource == "allcontacts-f" ? "All polar contacts, RS" :
                          "dataset???"
      if (clauses.length == 0) return datasetName

      return `${datasetName} with other filters`
    }

    function ensureLength(array: NumRange[], index: number) {
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
      if (["ranges", "basic"].includes(mode) && filter.sql.trim() == currentSqlQuery.trim()) {
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

    function getOrderByOptions(currentOption: string, allowed: string[] | null) {
      const virtualOpt = orderByOptions.find(x => x.id == currentOption) ? [] : [{ id: currentOption, expr: currentOption, label: _.truncate(currentOption, {length: 60}), title: "Custom sort expression " + currentOption }]
      if (allowed == null) {
        return [...orderByOptions, ...virtualOpt ]
      }
      else {
        return [...orderByOptions.filter(x => allowed.includes(x.id) || x.id == currentOption), ...virtualOpt ]
      }
    }

    async function setFr3dObservedBoundaries(sql: boolean) {
      const f = await filterfilter.defaultFilterLimits.value
      const newFilter = filterfilter.toNtFilter(f, metadata.pair_type.join("-"), null)
      const hbLengthLimits = metadata.atoms.map(([_, a, b, __]) =>
        a.includes("C") || b.includes("C") ? 4 : 3.8)
      newFilter.bond_length = hbLengthLimits.map(l => ({ min: null, max: l }))
      newFilter.datasource = filter.datasource
      newFilter.filtered = filter.filtered && !["fr3d-f", "allcontacts-f", "allcontacts-f"].includes(filter.datasource)
      if (sql) {
        newFilter.sql = makeSqlQuery(newFilter, getDataSourceTable(newFilter))
      }
      newFilter.filtered = filter.filtered
      newFilter.rotX = filter.rotX
      newFilter.orderBy = filter.orderBy
      newFilter.dna = filter.dna
      newFilter.resolution = filter.resolution
      filter = newFilter
    }

    let hasYawPitchRoll = false
    $: hasYawPitchRoll = Boolean(filter.yaw1 || filter.pitch1 || filter.roll1 || filter.yaw2 || filter.pitch2 || filter.roll2)

    function dataSourceChange(newDS: string) {
      if (newDS == "fr3d-f" || newDS == null) {
        filter = {...filter, datasource: undefined, filtered: true }
      }
      else if (newDS == "allcontacts-boundaries-f") {
        filter = {...filter, datasource: "allcontacts-f", filtered: true }
        setFr3dObservedBoundaries(false)
      } else {
        const filtered = newDS.endsWith('-f') || newDS.endsWith('-nf')
        filter = {...filter, datasource: newDS as any, filtered: filtered }
      }
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
          <input type="radio" checked={mode=="basic"} value="basic" name="editor_mode" on:change={modeChange}>
          Basic
        </label>
        <label class="radio" title="Filter by constraining the H-bond parameters.">
          <input type="radio" checked={mode=="ranges"} value="ranges" name="editor_mode" on:change={modeChange}>
          Parameter ranges
        </label>
        <label class="radio" title="Filter by anything you want using SQL.">
          <input type="radio" checked={mode=="sql"} value="sql" name="editor_mode" on:change={modeChange}>
          SQL
        </label>
    </div>

    {#if mode=="basic"}
    
    <div class="flex-columns" >
      <div class="column">
        <div class="field">
          <label class="label" for="ntfilter-data-source">Data source</label>
          <div class="control">
            <div class="select is-small">
              <select
                value={filter.datasource ?? 'fr3d-f'}
                id="ntfilter-data-source"
                on:change={ev => {
                  dataSourceChange(ev.currentTarget.value)
                }}
              >
                <option value="fr3d-f">FR3D, Representative Set</option>
                <option value="fr3d">FR3D, entire PDB</option>
                <option value="fr3d-nf">FR3D with nears, RS</option>
                <option value="fr3d-n">FR3D with nears, PDB</option>
                <option value="allcontacts-f">All polar contacts, RS</option>
                <option value="allcontacts-boundaries-f">All with boundaries, RS</option>
              </select>
            </div>
          </div>
        </div>
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

      </div>
      <div class="column">
        <div class="field">
          <label class="label" for="ntfilter-order-by">Order by</label>
          <div class="control">
            <div class="select is-small">
              <select bind:value={filter.orderBy} id="ntfilter-order-by">
                {#each getOrderByOptions(filter.orderBy, ["", "pdbidD", "rmsdA", "rmsdD"]) as opt}
                  <option value={opt.id} title={opt.title}>{opt.label}</option>
                {/each}
              </select>
            </div>
          </div>
        </div>

        <div class="control">
          <label class="checkbox" title="Rotate images 90° along X-axis to see the coplanarity">
            <input type="checkbox" checked={!!filter.rotX} on:change={e => filter = {...filter, rotX: e.currentTarget.checked }}>
            Rotate images
          </label>
        </div>
      </div>
      <div class="column" style="display: flex; flex-direction: column; justify-content: center">
        <button class="button is-warning" on:click={() => { filterBaseline = null; filter = defaultFilter() } }>Reset filters</button>

        {#if allowFilterBaseline && filterBaseline == null}
          {#if filter.datasource?.startsWith("allcontacts")}
            <button class="button" type="button" on:click={() => setBaseline({ ...defaultFilter(), datasource: filter.filtered ? "fr3d-f" : "fr3d" }, "ranges")}>
              Compare with FR3D
            </button>
          {:else}
            <!-- <button class="button" type="button" on:click={() => setBaseline({...filter}, mode)}>
              Compare with ???
            </button> -->
          {/if}
        {:else if filterBaseline != null}
          <button class="button is-warning" type="button" on:click={() => setBaseline(null, mode)}>
            Exit comparison
          </button>
        {/if}
      </div>
    </div>
    
    {:else if mode=="ranges"}
    <div class="flex-columns" >
        <div class="column">
          <div class="panel-title"></div>
          {#each bonds as bond, i}
            <div class="panel-field">{bond}</div>
          {/each}
          {#if hasYawPitchRoll}
            <div class="panel-title"></div>
            <div class="panel-field">Left to right</div>
            <div class="panel-field">Right to left</div>
          {/if}
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
            {#if hasYawPitchRoll}
              <h3 class="panel-title">Yaw</h3>
              <div class="panel-field field is-horizontal">
                <RangeEditor bind:range={filter.yaw1} step={1} min={-180} max={180} />
              </div>
              <div class="panel-field field is-horizontal">
                <RangeEditor bind:range={filter.yaw2} step={1} min={-180} max={180} />
              </div>
            {/if}

            {#if filter.coplanarity_angle}
              <h3 class="panel-title">Coplanarity angle</h3>
              <div class="panel-field field is-horizontal">
                <RangeEditor bind:range={filter.coplanarity_angle} step={1} min={-180} max={180} />
              </div>
            {/if}
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
            {#if hasYawPitchRoll}
              <h3 class="panel-title">Pitch</h3>
              <div class="panel-field field is-horizontal">
                <RangeEditor bind:range={filter.pitch1} step={1} min={-180} max={180} />
              </div>
              <div class="panel-field field is-horizontal">
                <RangeEditor bind:range={filter.pitch2} step={1} min={-180} max={180} />
              </div>
            {/if}

            {#if filter.coplanarity_edge_angle1}
              <h3 class="panel-title">Edge1/Plane2 angle</h3>
              <div class="panel-field field is-horizontal">
                <RangeEditor bind:range={filter.coplanarity_edge_angle1} step={1} min={-180} max={180} />
              </div>
            {/if}

            {#if filter.coplanarity_shift1}
              <h3 class="panel-title">Edge1/Plane2 distance</h3>
              <div class="panel-field field is-horizontal">
                <RangeEditor bind:range={filter.coplanarity_shift1} step={1} min={-180} max={180} />
              </div>
            {/if}
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
            {#if hasYawPitchRoll}
              <h3 class="panel-title">Roll</h3>
              <div class="panel-field field is-horizontal">
                <RangeEditor bind:range={filter.roll1} step={1} min={-180} max={180} />
              </div>
              <div class="panel-field field is-horizontal">
                <RangeEditor bind:range={filter.roll2} step={1} min={-180} max={180} />
              </div>
            {/if}

            {#if filter.coplanarity_edge_angle2}
              <h3 class="panel-title">Coplanarity Edge2/Plane1</h3>
              <div class="panel-field field is-horizontal">
                <RangeEditor bind:range={filter.coplanarity_edge_angle2} step={1} min={-180} max={180} />
              </div>
            {/if}
            {#if filter.coplanarity_shift2}
              <h3 class="panel-title">Edge2/Plane1 distance</h3>
              <div class="panel-field field is-horizontal">
                <RangeEditor bind:range={filter.coplanarity_shift2} step={1} min={-180} max={180} />
              </div>
            {/if}
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

          <div class="field">
            <label class="label" for="ntfilter-data-source">Data source</label>
            <div class="control">
              <div class="select is-small">
                <select
                  value={filter.datasource ?? 'fr3d-f'}
                  id="ntfilter-data-source"
                  on:change={ev => {
                    dataSourceChange(ev.currentTarget.value)
                  }}
                >
                  <option value="fr3d-f">FR3D, Representative Set</option>
                  <option value="fr3d">FR3D, entire PDB</option>
                  <option value="fr3d-nf">FR3D with nears, RS</option>
                  <option value="fr3d-n">FR3D with nears, PDB</option>
                  <option value="allcontacts-f">All polar contacts, RS</option>
                  <option value="allcontacts-boundaries-f">All with boundaries, RS</option>
                </select>
              </div>
            </div>
          </div>

          {#if [].includes(filter.datasource)}
          <div class="control">
            <label class="checkbox" title="Filter out redundant nucleotides or nucleoties with bad something TODO">
              <input type="checkbox" checked={filter.filtered} on:change={e => filter = {...filter, filtered: e.currentTarget.checked }}>
              Representative set only
            </label>
          </div>
          {/if}

          <div class="field has-addons">
            {#if filter.resolution?.min != null}
              <div class="control">
                <input class="input is-small num-input" style="max-width:4rem" type="number" step="0.1" min=0 max={filter.filtered ? 3.5 : 20} placeholder="Min" value={filter.resolution?.min ?? 0} on:change={ev => { filter.resolution ??= {}; filter.resolution.min = tryParseNum(ev.currentTarget.value)} }>
              </div>
            {/if}
            <label class="label" for="ntfilter-resolution">{#if filter.resolution?.min != null}&nbsp;≤ {/if}Resolution ≤&nbsp;</label>
            <div class="control">
              <input class="input is-small num-input" style="max-width:4rem" type="number" step="0.1" min=0 max={filter.filtered ? 3.5 : 20} placeholder={filter.filtered ? '3.5' : ''} value={filter.resolution?.max ?? ''} on:change={ev => { filter.resolution ??= {}; filter.resolution.max = tryParseNum(ev.currentTarget.value)}}>
            </div>
            &nbsp;Å
          </div>
          {#if filter.filtered && filter.resolution?.max && filter.resolution?.max > 3.5}
            <p class="help is-danger">Representative set only<br> contains structures ≤3.5 Å</p>
          {/if}
        </div>

        <div class="column">
          <div class="field">
            <label class="label" for="ntfilter-order-by">Order by</label>
            <div class="control">
              <div class="select is-small">
                <select bind:value={filter.orderBy} id="ntfilter-order-by">
                  {#each getOrderByOptions(filter.orderBy, null) as opt}
                    <option value={opt.id} title={opt.title}>{opt.label}</option>
                  {/each}
                </select>
              </div>
            </div>
          </div>

          {#if allowFilterBaseline}
          <div class="field">
            <label class="label" for="ntfilter-order-by">Comparison baseline</label>
            <div class="control">
              <div class="buttons has-addons">
                {#if filterBaseline == null}
                <button class="button is-small" on:click={() => setBaseline({ ...defaultFilter(), datasource: filter.filtered ? "fr3d-f" : "fr3d" }, "ranges")}
                  title="Sets the current filters as the filter baseline, allowing you to change some parameters and observe the changed">
                  FR3D
                </button>
                <button class="button is-small" on:click={() => setBaseline({...filter}, mode)}
                  title="Compares the current basepair selection with that determined by FR3D">
                  Set to this
                </button>
                {:else}
                  <button class="button is-warning" type="button" on:click={() => setBaseline(null, mode)}
                    title="Exits comparison mode, removed the baseline">
                    ❌ Reset
                  </button>
                  <!-- <button class="button" type="button" on:click={() => setBaseline(filter, mode)}>
                    Set to this
                  </button> -->

                {/if}
              </div>
            </div>
          </div>
          {/if}
          
          <div class="control">
            <label class="checkbox" title="Rotate images 90° along X-axis to see the coplanarity">
              <input type="checkbox" checked={!!filter.rotX} on:change={e => filter = {...filter, rotX: e.currentTarget.checked }}>
              Rotate images
            </label>
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
      <div style="float: right; margin-right: 2rem">
        <div class="control" style="display: inline">
          <label class="checkbox" title="Rotate images 90° along X-axis to see the coplanarity">
            <input type="checkbox" checked={!!filter.rotX} on:change={e => filter = {...filter, rotX: e.currentTarget.checked }}>
            Rotate images
          </label>
        </div>
      </div>
      {#if allowFilterBaseline}
        <div style="display: flex; align-items: center; gap: 1rem">
        {#if filter.datasource?.startsWith("allcontacts") && metadata != null}
          <button class="button" type="button" on:click={()=> {
            setFr3dObservedBoundaries(true)
          } }>
            Constrain to FR3D observed ranges
          </button>
        {/if}
        {#if filterBaseline == null}
          <button class="button" type="button" on:click={() => setBaseline({...filter}, mode)}>
            Set as baseline
          </button>
          <button class="button" type="button" on:click={() => setBaseline({ ...defaultFilter(), datasource: filter.filtered ? "fr3d-f" : "fr3d" }, "ranges")}>
            Compare with FR3D
          </button>
        {:else}
          <button class="button is-warning" type="button" on:click={() => setBaseline(null, mode)}>
            Exit comparison
          </button>
          <button class="button" type="button" on:click={() => setBaseline(filter, mode)}>
            Set current query as baseline
          </button>
          {#if filterBaseline.sql == null}
          <button class="button" type="button" on:click={() => filterBaseline = {...filterBaseline, sql: makeSqlQuery(filterBaseline, getDataSourceTable(filterBaseline)) }}>
            Edit baseline
          </button>
          {/if}
          {#if comparisonMode != null}
          <div class="select">
            <select bind:value={comparisonMode}>
              <option value="union">Show all matches</option>
              <option value="difference">Only differences</option>
              <option value="new">Absent in baseline</option>
              <option value="missing">Only in baseline</option>
            </select>
          </div>
          {/if}
          <span>Comparing with {formatFilter(filterBaseline, filter)}</span>
        {/if}
        </div>
      {/if}

      {#if filterBaseline != null && filterBaseline.sql != null}
        <div class="field">
          <label class="label is-medium" for="filter-baseline-sql-textarea">Comparing with baseline query:</label>
          <div class="control">
            <textarea class="textarea sql-editor" id="filter-baseline-sql-textarea" bind:value={filterBaseline.sql} style="width: 100%; font-color: #bbbbbb"></textarea>
          </div>
        </div>
      {/if}
    {/if}
</div>
