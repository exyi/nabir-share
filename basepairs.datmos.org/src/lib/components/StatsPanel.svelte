<script lang="ts">
    import { fillStatsLegends, getColumnLabel, hideColumn, type HistogramSettingsModel, type KDE2DSettingsModel, type StatisticsSettingsModel } from "$lib/dbModels";
    import type * as arrow from 'apache-arrow'
    import HistogramPlot from "./HistogramPlot.svelte";
    import Density2DPlot from "./Density2DPlot.svelte";
    import { variance } from "d3";

    export let settings: StatisticsSettingsModel | null = null
    export let data: arrow.Table | null = null
    export let availableSchema: arrow.Schema | null = null
    export let metadata: any = null

    export let allowEditing: boolean = true
    $: {
      console.log("stat panels update: ", settings, settings.panels[0]?.variables, settings.panels[0]?.variables[0], settings.panels[0]?.variables[0]?.label)
    }

    let modalDialog: HTMLDialogElement
    let rootElement: HTMLElement
    let optionsIndex: number | null = null
    
    let options: HistogramSettingsModel | KDE2DSettingsModel | undefined
    $: options = settings?.panels[optionsIndex]
    let optionsAllowFilters = true
    $: optionsAllowFilters = !!options?.variables?.some(v => !!v.filterSql)
    let availableColumns
    $: availableColumns = data && metadata ? (availableSchema ?? data.schema).fields.filter(f => !hideColumn(f.name, metadata)).map(f => {
      const [label, tooltip] = getColumnLabel(f.name, metadata) ?? [f.name, null]
      return { column: f.name, label, tooltip: `${f.name}: ${f.type}${tooltip ? ' - ' + tooltip : ''}`, type: ""+f.type, isNumber: String(f.type).startsWith("Float") /* || String(f.type).startsWith("Int") */ }
    }) : []

    function editPanel(panelIndex: number) {
      modalDialog.showModal()
      const sender = rootElement.querySelector(`div.stat-panel[data-statpanel-index="${panelIndex}"] .options-icon`)
      optionsIndex = panelIndex
      const pos = sender.getBoundingClientRect()
      modalDialog.style.top = (window.scrollY+pos.top) + "px"
      modalDialog.style.left = `calc(min(${window.scrollX + window.innerWidth}px - var(--stat-panel-width), ${window.scrollX+pos.left}px))`
    }

</script>

<style>
    .stat-panels {
        display: flex;
        flex-direction: row;
        flex-wrap: wrap;
        justify-content: space-evenly;
    }
    .stat-panel {
        display: flex;
        flex-direction: column;
        /* align-items: center; */
        margin: 0.5rem;
    }
    .stat-panel > * {
        width: 100%;
        width: var(--stat-panel-width, 33vw - 2rem);
    }
    @media (max-width: 800px) {
        .stat-panels {
          --stat-panel-width: calc(100vw - 2rem);
        }
    }
    @media (max-width: 1200px) and (min-width: 800px) {
        .stat-panels {
          --stat-panel-width: calc(50vw - 2rem);
        }
    }
    @media (min-width: 1200px) {
        .stat-panels {
          --stat-panel-width: calc(33vw - 2rem);
        }
    }
    /* @media (min-width: 2000px) {
        .stat-panel {
            width: calc(25vw - 2rem);
        }
    } */
    .title {
        text-align: center;
        /* margin-bottom: -10px; */
    }

    .options-icon {
        float:right;
        font-size: 1rem;
        margin-right: 0.5rem;
        cursor: pointer;
        opacity: 0.5;
        background-color: #cccccc;
        border-radius: 0.5rem;
        padding: 0.25rem;
        border: none;
    }
    .options-icon:hover {
        opacity: 1;
    }

    .add-new-panel-button-hack {
      overflow: visible;
      /* height: 0px; */
      width: 100%;
      margin-top: -20px;
      opacity: 0.5;
      text-align: right;
    }
    .add-new-panel-button-hack:hover {
      opacity: 1;
    }

    dialog {
      position: absolute;
      width: var(--stat-panel-width, 33vw - 2rem);
      margin: 0;
    }

</style>

<div class="stat-panels" bind:this={rootElement}>
  {#each fillStatsLegends(settings, metadata).panels as panel, panelIndex}
    <div class="stat-panel" data-statpanel-index={panelIndex}>
      <h4 class="title is-6">
        {panel.title || '.'}
        {#if allowEditing}
        <button class="options-icon" on:click={(ev) => {
          editPanel(panelIndex)
        }}>⚙️</button>
        {/if}
      </h4>
      {#if panel.type == "histogram"}
        <HistogramPlot data={data} settings={panel} />
      {:else if panel.type == 'kde2d'}
        <Density2DPlot data={data} settings={panel} />
      {:else}
        Unknown panel type: {panel.type}
      {/if}
      {#if allowEditing && panelIndex == settings.panels.length - 1}
      <div class="add-new-panel-button-hack">
        <button class="button is-success" title="Add new stat panel" on:click={() => {
          settings.panels = [
            ...settings.panels,
            { title: "New plot", type: "histogram", variables: [{ column: "hb_0_length", label: "" }] }
          ]
          setTimeout(() => editPanel(settings.panels.length - 1), 200)
        }}>New</button>
      </div>
      {/if}
    </div>
  {/each}
  {#if allowEditing && settings.panels.length == 0}
  <div class="stat-panel" style="text-align: center">
    <button class="button is-success" title="Add new stat panel" on:click={() => {
      settings.panels = [
        ...settings.panels,
        { title: "New plot", type: "histogram", variables: [{ column: "hb_0_length", label: "" }] }
      ]
      setTimeout(() => editPanel(settings.panels.length - 1), 200)
    }}>Create new plot</button>
  </div>
  {/if}
  <dialog class="panel-options" bind:this={modalDialog} on:click={ev => {
    const dialog = ev.currentTarget;
    const rect = dialog.getBoundingClientRect();
    const isInDialog = (rect.top <= ev.clientY && ev.clientY <= rect.top + rect.height &&
      rect.left <= ev.clientX && ev.clientX <= rect.left + rect.width);
    if (!isInDialog && !["SELECT", "OPTION"].includes(ev.target.tagName)) {
      dialog.close();
    }
  }}>
    <form>
    {#if optionsIndex != null}
    <h5 class="title is-6">Options for panel {optionsIndex} {options.type}</h5>
    <div class="field has-addons">
      <div class="control">
        <div class="select">
          <select
            value={options.type}
            on:change={ev => {
              const newVal = ev.currentTarget.value
              if (newVal == options.type) {
                return
              }
              else if (newVal == "histogram") {
                settings.panels[optionsIndex].type = "histogram"
              } else if (newVal == "kde2d") {
                settings.panels[optionsIndex].type = "kde2d"
                while (options.variables.length < 2) {
                  options.variables.push({ column: "", label: "" })
                }
              }
            }}>
            <option value="histogram">📊 Histogram</option>
            <option value="kde2d">🌌 2D Scatterplot</option>
          </select>
        </div>
      </div>
      <div class="control is-expanded">
        <input class="input" type="text" style="font-weight: bold;" bind:value={settings.panels[optionsIndex].title}>
      </div>
    </div>
    {#if settings.panels[optionsIndex].type == "histogram"}
      {#each settings.panels[optionsIndex].variables as v}
        <div class="field has-addons">
          <div class="control">
            <div class="select">
              <select bind:value={v.column} on:change={e => v.label = ""}>
                {#each availableColumns.filter(c => c.isNumber) as col}
                  <option value={col.column} title="{col.tooltip}">{col.label}</option>
                {/each}
              </select>
            </div>
          </div>
          {#if optionsAllowFilters}
          <div class="control">
            <div class="input">WHERE</div>
          </div>
          <div class="control is-expanded">
            <input class="input" type="text" bind:value={v.filterSql} placeholder="">
          </div>
          {/if}
          <div class="control">
            <button class="button is-danger" title="Remove variable" on:click={() => {
              settings.panels[optionsIndex].variables = settings.panels[optionsIndex].variables.filter(x => x != v)
            }}>⨯</button>
          </div>
        </div>
      {/each}
      <div class="field is-grouped is-grouped-right">
        <button class="button is-success" title="Add new variable" on:click={() => {
          settings.panels[optionsIndex].variables = [
            ...settings.panels[optionsIndex].variables,
            { column: "", label: "" }
          ]
        }}>＋</button>
      </div>
    {:else if settings.panels[optionsIndex].type == "kde2d"}
      
      {#each settings.panels[optionsIndex].variables.slice(0, 2) as v, vi}
      <label class="label" for="input_statpanel_variable_{vi}">{["X-axis (horizontal)", "Y-axis (vertical)"][vi]}</label>
      <div class="field has-addons">
        <div class="control">
          <div class="select">
            <select bind:value={v.column} id="input_statpanel_variable_{vi}" on:change={e => v.label = ""}>
              {#each availableColumns.filter(c => c.isNumber) as col}
                <option value={col.column} title="{col.tooltip}">{col.label}</option>
              {/each}
            </select>
          </div>
        </div>
        {#if optionsAllowFilters}
        <div class="control">
          <div class="input">WHERE</div>
        </div>
        <div class="control is-expanded">
          <input class="input" type="text" bind:value={v.filterSql} placeholder="">
        </div>
        {/if}
      </div>
      {/each}

    {/if}

    <div class="field is-grouped">
      <button type="submit" class="control button is-primary" on:click={ev => {
        modalDialog.close()
        optionsIndex = null
        ev.stopPropagation()
      }}>Ok, close</button>
      <button class="control button is-danger" on:click={ev => {
        settings.panels = settings.panels.filter((x, i) => i != optionsIndex)
        modalDialog.close()
        optionsIndex = null
        ev.stopPropagation()
      }}>Delete</button>

      <button class="control button is-light" disabled={optionsIndex == 0} on:click={ev => {
        if (optionsIndex > 0) {
          const next = optionsIndex - 1
          const x = settings.panels[optionsIndex]
          settings.panels[optionsIndex] = {...settings.panels[next]}
          settings.panels[next] = {...x}
          modalDialog.close()
          optionsIndex = null
          editPanel(next)
        }
        ev.stopPropagation()
      }}>↩</button>
      <button class="control button is-light" disabled={optionsIndex >= settings.panels.length - 1} on:click={ev => {
        if (optionsIndex < settings.panels.length - 1) {
          const next = optionsIndex + 1
          const x = settings.panels[optionsIndex]
          settings.panels[optionsIndex] = {... settings.panels[next]}
          settings.panels[next] = {...x}
          modalDialog.close()
          optionsIndex = null
          editPanel(next)
        }
        ev.stopPropagation()
      }}>↪</button>
    </div>
    {/if}
    </form>
  </dialog>
</div>
