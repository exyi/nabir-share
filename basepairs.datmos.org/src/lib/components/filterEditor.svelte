<script lang="ts">
	import { filterToSqlCondition, makeSqlQuery, type NucleotideFilterModel } from "$lib/dbLayer";

    export let filter: NucleotideFilterModel
    export let selectingFromTable: string | null = null
    export let mode: "ranges" | "sql" = "ranges"

    let bonds = ["0", "1", "2"]
    let ranges = {
        length: bonds.map(_ => ({ min: "", max: "" })),
        accAngle: bonds.map(_ => ({ min: "", max: "" })),
        donAngle: bonds.map(_ => ({ min: "", max: "" })),
    }
    function updateFilter() {
        function convertRange(x: any) {
            return { min: (isNaN(Number.parseFloat(x.min)) ? null : Number.parseFloat(x.min)),
                     max: (isNaN(Number.parseFloat(x.max)) ? null : Number.parseFloat(x.max)) }
        }
        filter = { ...filter,
            bond_length: bonds.map((_, i) => convertRange(ranges.length[i])),
            bond_acceptor_angle: bonds.map((_, i) => convertRange(ranges.accAngle[i])),
            bond_donor_angle: bonds.map((_, i) => convertRange(ranges.donAngle[i]))
        }
    }

    $:{
        ranges
        updateFilter()
    }

    function modeChange(e: Event & { currentTarget:EventTarget & HTMLInputElement }) {
      mode = e.currentTarget.value as any

      if (mode=="sql" && !filter.sql) {
        filter = {...filter, sql: makeSqlQuery(filter, selectingFromTable) }
      }
      if (mode == "ranges" && filter.sql.trim() == makeSqlQuery(filter, selectingFromTable).trim()) {
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
  }
</style>

<div>
    <div class="control">
        <label class="radio">
          <input type="radio" checked={mode=="ranges"} value="ranges" name="editor_mode" on:change={modeChange}>
          Ranges
        </label>
        <label class="radio">
          <input type="radio" checked={mode=="sql"} value="sql" name="editor_mode" on:change={modeChange}>
          SQL
        </label>
    </div>
    
    {#if mode=="ranges"}
    <div class="columns" >
        <div class="column">
            <h3 class="panel-title">Length</h3>
            {#each bonds as bond, i}
            <div class="field is-horizontal">
                <div class="field-label">Bond {bond}</div>
                <div class="field-body">
                  <div class="field is-expanded">
                    <div class="field has-addons">
                      <p class="control">
                        <input class="input is-small" type="number" step="0.1" min=0 max=6 placeholder="Min" bind:value={ranges.length[i].min}>
                      </p>
                      <p class="control is-expanded">
                        <input class="input is-small" type="number" step="0.1" min=0 max=6 placeholder="Max" bind:value={ranges.length[i].max}>
                      </p>
                    </div>
                  </div>
                </div>
              </div>
            {/each}
        </div>
        <div class="column">
            <h3 class="panel-title">Donor Angle</h3>
            {#each bonds as bond, i}
            <div class="field is-horizontal">
                <div class="field-label">Bond {bond}</div>
                <div class="field-body">
                  <div class="field is-expanded">
                    <div class="field has-addons">
                      <p class="control">
                        <input class="input is-small" type="number" step="5" min=0 max=360 placeholder="Min" bind:value={ranges.donAngle[i].min}>
                      </p>
                      <p class="control is-expanded">
                        <input class="input is-small" type="number" step="5" min=0 max=360 placeholder="Max" bind:value={ranges.donAngle[i].max}>
                      </p>
                    </div>
                  </div>
                </div>
              </div>
            {/each}
        </div>

        <div class="column">
            <h3 class="panel-title">Acceptor Angle</h3>
            {#each bonds as bond, i}
            <div class="field is-horizontal">
                <div class="field-label">Bond {bond}</div>
                <div class="field-body">
                  <div class="field is-expanded">
                    <div class="field has-addons">
                      <p class="control">
                        <input class="input is-small" type="number" step="5" min=0 max=360 placeholder="Min" bind:value={ranges.accAngle[i].min}>
                      </p>
                      <p class="control is-expanded">
                        <input class="input is-small" type="number" step="5" min=0 max=360 placeholder="Max" bind:value={ranges.accAngle[i].max}>
                      </p>
                    </div>
                  </div>
                </div>
              </div>
            {/each}
        </div>

        <div class="column">
          <div class="control">
            <label class="radio">
              <input type="radio" name="rna_dna_mode" value="rna" checked={filter.dna == false} on:change={dnaRnaChange}>
              RNA
            </label>
            <label class="radio">
              <input type="radio" name="rna_dna_mode" value="dna" checked={filter.dna == true} on:change={dnaRnaChange}>
              DNA
            </label>
            <label class="radio">
              <input type="radio" name="rna_dna_mode" value="both" checked={filter.dna == null} on:change={dnaRnaChange}>
              Both
            </label>
          </div>

          <div class="field">
            <label class="label">Order by</label>
            <div class="control">
              <div class="select is-small">
                <select bind:value={filter.orderBy}>
                  <option value="">pdbid</option>
                  <option value="pdbid desc, model desc, chain1 desc, nr1 desc">pdbid descending</option>
                  <option value="resolution, pdbid, model, chain1, nr1">resolution</option>
                </select>
              </div>
            </div>
          </div>
        </div>

    </div>

    {:else if mode=="sql"}
      <div>
        <textarea class="textarea" bind:value={filter.sql} style="width: 100%;"></textarea>
      </div>
    {/if}
</div>
