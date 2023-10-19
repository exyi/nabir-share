<script lang="ts">
	import type { NucleotideFilterModel } from "$lib/dbLayer";

    export let filter: NucleotideFilterModel
    export let isPlainSql: boolean

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
</script>

<div>
    <div class="control">
        <label class="radio">
          <input type="radio" name="editor_mode">
          Ranges
        </label>
        <label class="radio">
          <input type="radio" name="editor_mode">
          SQL
        </label>
        <!-- <label class="radio">
          <input type="radio" name="editor_mode" disabled>
          Maybe
        </label> -->
    </div>

    <div class="columns">
        <div class="column">
            <h3>length</h3>
            {#each bonds as bond, i}
            <div class="field is-horizontal">
                <div class="field-label">Bond {bond}</div>
                <div class="field-body">
                  <div class="field is-expanded">
                    <div class="field has-addons">
                      <p class="control">
                        <input class="input" type="number" step="0.1" min=0 max=6 placeholder="Min" bind:value={ranges.length[i].min}>
                      </p>
                      <p class="control is-expanded">
                        <input class="input" type="number" step="0.1" min=0 max=6 placeholder="Max" bind:value={ranges.length[i].max}>
                      </p>
                    </div>
                  </div>
                </div>
              </div>
            {/each}
        </div>
        <div class="column">
            <h3>donor angle</h3>
            {#each bonds as bond, i}
            <div class="field is-horizontal">
                <div class="field-label">Bond {bond}</div>
                <div class="field-body">
                  <div class="field is-expanded">
                    <div class="field has-addons">
                      <p class="control">
                        <input class="input" type="number" step="0.1" min=0 max=6 placeholder="Min" bind:value={ranges.length[i].min}>
                      </p>
                      <p class="control is-expanded">
                        <input class="input" type="number" step="0.1" min=0 max=6 placeholder="Max" bind:value={ranges.length[i].max}>
                      </p>
                    </div>
                  </div>
                </div>
              </div>
            {/each}
        </div>

        <div class="column">
            <h3>acceptor angle</h3>
            {#each bonds as bond, i}
            <div class="field is-horizontal">
                <div class="field-label">Bond {bond}</div>
                <div class="field-body">
                  <div class="field is-expanded">
                    <div class="field has-addons">
                      <p class="control">
                        <input class="input" type="number" step="0.1" min=0 max=6 placeholder="Min" bind:value={ranges.length[i].min}>
                      </p>
                      <p class="control is-expanded">
                        <input class="input" type="number" step="0.1" min=0 max=6 placeholder="Max" bind:value={ranges.length[i].max}>
                      </p>
                    </div>
                  </div>
                </div>
              </div>
            {/each}
        </div>

    </div>
</div>
