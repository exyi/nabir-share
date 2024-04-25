<script lang="ts">
    import type { HistogramSettingsModel } from "$lib/dbModels";
    import * as d3 from "d3";
    import type * as arrow from 'apache-arrow'
	import * as a from "$lib/arrowUtils";
	import { onMount } from "svelte";
	import { bind } from "svelte/internal";
    export let settings: HistogramSettingsModel | null = null
    export let data: arrow.Table | null = null

    let svg: SVGGElement
    const widthPx = 640
    const heightPx = 480
    const marginPx = {top: 10, right: 10, bottom: 20, left: 40}
    const colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999']

    let hiddenColumns: number[] = []

    // function createRoundedAxis(min, max, buckets = 20) {
    //     const ticks = d3.scaleLinear().domain([min, max]).ticks(buckets)
    //     const tickDiff = ticks[2] - ticks[1]
    //     const min2 = ticks[1] - tickDiff
    //     const max2 = ticks[ticks.length - 2] + tickDiff
    //     return d3.scaleLinear().domain([min2, max2]).ni
    // }

    function makeD3Bins(data: (a.ColumnHelper | null)[], min: number, max: number, binCount: number) {
        const binWidth = (max - min) / binCount
        const binData = data.map(c => c == null ? [] : Array.from(a.binArrays(c.data, min, max, binCount)))
        const d3Bins = binData.map(b => b.map((count, i) => ({
            count,
            x0: min + i * binWidth,
            x1: min + (i + 1) * binWidth
        })))
        return { d3Bins, max: Math.max(...binData.flatMap(b => Math.max(...b))) }
    }

    function rerender() {
        if (!settings) return
        const startTime = performance.now()

        const variables = settings.variables //.filter((_, i) => !hiddenColumns.includes(i))
        const columns = variables.map((v, i) => hiddenColumns.includes(i) || !v.column ? null : a.tryGetColumnHelper(data, v))
        let min = Math.min(...columns.filter(c => c != null).map(c => c.data.reduce((a, b) => Math.min(a, b), 1e50)))
        let max = Math.max(...columns.filter(c => c != null).map(c => c.data.reduce((a, b) => Math.max(a, b), -1e50)))
        let defaultBinWidth: null | number = null
        let defaultBins: null | number = null
        if (settings.variables[0]?.column.endsWith("_length")) {
            min = Math.min(2.5, min)
            max = Math.max(3.5, max)
            defaultBinWidth = 0.05
        }

        const xAxis = d3.scaleLinear().domain([min, max])
            .nice(settings.bins ?? defaultBins ?? 40)
            .range([marginPx.left, widthPx - marginPx.right])
        defaultBins ??= defaultBinWidth != null ? Math.ceil((xAxis.domain()[1] - xAxis.domain()[0]) / defaultBinWidth) : 40
        const binCount = settings.bins ?? defaultBins


        // const bins = d3.bin().thresholds(xAxis.ticks(settings?.bins ?? 40)).domain(xAxis.domain() as [number, number])
        // const binData = columns.map(c => bins(c.data))
        const bins = makeD3Bins(columns, xAxis.domain()[0], xAxis.domain()[1], binCount)
        const yAxis = d3.scaleLinear()
            .domain([0, Math.max(1, bins.max)])
            .range([heightPx - marginPx.bottom, marginPx.top])

        const node = d3.select(svg)
        node.selectAll("*").remove()
        node.append("g")
            .attr("transform", `translate(0,${heightPx - marginPx.bottom})`)
            .call(d3.axisBottom(xAxis))
        node.append("g")
            .attr("transform", `translate(${marginPx.left},0)`)
            .call(d3.axisLeft(yAxis))
        for (let i = 0; i < columns.length; i++) {
            if (hiddenColumns.includes(i))
                continue
            const c = columns[i]
            const b = bins.d3Bins[i]
            const color = colors[i] ?? '#000000'
            node.append("g")
                .selectAll("rect")
                .data(b)
                .enter()
                .append("rect")
                .attr("transform", function(d) { return `translate(${xAxis(d.x0)},${yAxis(d.count)})`; })
                .attr("width", function(d) { return xAxis(d.x1) - xAxis(d.x0) - 1; })
                .attr("height", function(d) { return heightPx - yAxis(d.count) - marginPx.bottom; })
                .style("fill", d3.color(color).brighter(0.5).hex())
                .style("stroke", d3.color(color).darker(0.5).hex())
                .style("fill-opacity", 0.6)
                .append("title")
                .text(d => `${d.x0} <= [${variables[i]?.label}] <= ${d.x1}: ${d.count}`)
        }
        // for (let colI = 0; colI < columns.length; colI++) {
        //     node.append("circle").attr("cx",300).attr("cy", ).attr("r", 6).style("fill", "#69b3a2")
        //     node.append("circle").attr("cx",300).attr("cy",60).attr("r", 6).style("fill", "#404080")
        //     svg.append("text").attr("x", 320).attr("y", 30).text("variable A").style("font-size", "15px").attr("alignment-baseline","middle")
        //     svg.append("text").attr("x", 320).attr("y", 60).text("variable B").style("font-size", "15px").attr("alignment-baseline","middle")
        // }
        console.log(`histogram rendered in ${performance.now() - startTime}ms`)
    }

    onMount(() => {
        rerender()
    })

    $: {
        if (settings && data && svg && hiddenColumns) {
            rerender()
        }
    }
</script>

<style>
    .histogram-plot {
        position: relative;
        /* width: 100%; */
    }
    svg {
        width: 100%;
        height: 100%;
    }
</style>

<div class="histogram-plot">
    <svg viewBox="0 0 {widthPx} {heightPx}">
        <g bind:this={svg}></g>
        <g>
        {#each settings?.variables ?? [] as v, i}
        {#if v.label != null}
            {@const maxLblLength = Math.max(...(settings?.variables ?? []).map(v => v.label?.length ?? 0))}

            <g opacity={hiddenColumns.includes(i) ? 0.5 : 1} on:click={() => {
                if (hiddenColumns.includes(i)) {
                    hiddenColumns = hiddenColumns.filter(x => x!=i)
                } else {
                    hiddenColumns = [ ...hiddenColumns, i ]
                }
            }} cursor="pointer">
                <circle cx={widthPx - marginPx.right - maxLblLength * 8 - 20} cy={marginPx.top + 30 * i + 6} r="6" style="fill: {d3.color(colors[i]).brighter(0.5).hex() ?? '#000000'}; fill-opacity: 0.6; stroke: {d3.color(colors[i]).darker(0.5).hex()}"/>
                <text x={widthPx - marginPx.right - maxLblLength * 8 - 10} y={marginPx.top + 30 * i + 12} text-anchor="start" alignment-baseline="middle">
                    <title>{v.column} {#if v.filterSql}WHERE {v.filterSql}{/if}</title>
                    {v.label}
                </text>
                <!-- <text x={widthPx - marginPx.right} y={marginPx.top + 30 * i + 12} text-anchor="end" alignment-baseline="middle" cursor="pointer">
                    {#if hiddenColumns.includes(i)}🐵{:else}🙈{/if}
                </text> -->
            </g>
        {/if}
        {/each}
        </g>
    </svg>
</div>
