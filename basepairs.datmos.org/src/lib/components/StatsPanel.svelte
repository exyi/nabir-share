<script lang="ts">
    import type { HistogramSettingsModel, StatisticsSettingsModel } from "$lib/dbModels";
    import type * as arrow from 'apache-arrow'
	import HistogramPlot from "./HistogramPlot.svelte";
	import Density2DPlot from "./Density2DPlot.svelte";

    export let settings: StatisticsSettingsModel | null = null
    export let data: arrow.Table<any> | null = null

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
        align-items: center;
        margin: 0.5rem;
    }
    @media (max-width: 800px) {
        .stat-panel {
            width: calc(100vw - 2rem);
        }
    }
    @media (max-width: 1200px) and (min-width: 800px) {
        .stat-panel {
            width: calc(50vw - 2rem);
        }
    }
    @media (min-width: 1200px) {
        .stat-panel {
            width: calc(33vw - 2rem);
        }
    }
</style>

<div class="stat-panels">
    {#each settings.panels as panel, panelIndex}
    <div class="stat-panel">
      <h4>{panel.title || '.'}</h4>
      {#if panel.type == "histogram"}
        <HistogramPlot data={data} settings={panel} />
      {:else if panel.type == 'kde2d'}
        <Density2DPlot data={data} settings={panel} />
      {:else}
        Unknown panel type: {panel.type}
      {/if}
    </div>
  {/each}
</div>
