<script lang="ts">
    import type { KDE2DSettingsModel } from "$lib/dbModels";
    import * as d3 from "d3";
    import type * as arrow from 'apache-arrow'
	import * as a from "$lib/arrowUtils";
	import { onMount } from "svelte";
	import { bind } from "svelte/internal";
    export let settings: KDE2DSettingsModel | null = null
    export let data: arrow.Table<any> | null = null

    let svg: SVGGElement
    const widthPx = 640
    const heightPx = 480
    const marginPx = {top: 10, right: 10, bottom: 32, left: 40}
    const colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999']

    function makeD3Bins(data: a.ColumnHelper[], min: number, max: number, binCount: number) {
        const binWidth = (max - min) / binCount
        const binData = data.map(c => a.binArrays(c.data, min, max, binCount))
        const d3Bins = binData.map(b => Array.from(b).map((count, i) => ({
            count,
            x0: min + i * binWidth,
            x1: min + (i + 1) * binWidth
        })))
        return { d3Bins, max: Math.max(...binData.flatMap(b => Math.max(...b))) }
    }

    function getColor(column: string) {
        let m
        if (m = /^hb_(\d+)/.exec(column)) {
            return colors[parseInt(m[1]) % colors.length]
        } else {
            return null
        }
    }

    function rerender() {
        const startTime = performance.now()
        if (!settings) return
        if (settings.variables.length != 2) {
            console.error("KDE2DPlot requires exactly two variables")
            return
        }

        const color = getColor(settings.variables[0].column) ?? getColor(settings.variables[1].column) ?? '#0000ff'

        const columns = settings.variables.map(v => a.getColumnHelper(data, v, settings.variables)) as [a.ColumnHelper, a.ColumnHelper]
        const x = columns[0].data,
              y = columns[1].data
        if (x.length != y.length) {
            console.error("KDE2DPlot requires variables to have the same length", settings.variables, columns)
            return
        }
        if (x.length <= 5) {
            return
        }

        const xAxis = d3.scaleLinear()
            .domain([
                x.reduce((a, b) => Math.min(a, b), 1e50),
                x.reduce((a, b) => Math.max(a, b), -1e50)
            ])
            .nice(40)
            .range([marginPx.left, widthPx - marginPx.right])
        const yAxis = d3.scaleLinear()
            .domain([
                y.reduce((a, b) => Math.min(a, b), 1e50),
                y.reduce((a, b) => Math.max(a, b), -1e50)
            ])
            .nice(40)
            .range([heightPx - marginPx.bottom, marginPx.top])

        const xScaled = x.map(x => xAxis(x))
        const yScaled = y.map(y => yAxis(y))
        const xStd = a.stddev(x)
        const yStd = a.stddev(y)

        const bandwidth = 0.5 * ((x.length)**(-1./(2+4))) * Math.max(xStd * Math.abs(xAxis.range()[1] - xAxis.range()[0]), yStd * Math.abs(yAxis.range()[1] - yAxis.range()[0]))

        const contours = d3.contourDensity<number>()
            .x(i => xScaled[i])
            .y(i => yScaled[i])
            .size([widthPx - marginPx.left - marginPx.right, heightPx - marginPx.top - marginPx.bottom])
            .bandwidth(bandwidth)
            .thresholds(40)
            (a.arange(0, xScaled.length) as any)
            
        const node = d3.select(svg)
        node.selectAll("*").remove()
        node.append("g")
            .attr("transform", `translate(0,${heightPx - marginPx.bottom})`)
            .call(d3.axisBottom(xAxis))
        node.append("g")
            .attr("transform", `translate(${marginPx.left},0)`)
            .call(d3.axisLeft(yAxis))
            .append("text")
            .attr("x", "2")
            .attr("y", "18")
            .attr("dy", ".71em")
            .text(settings.variables[1].label || settings.variables[1].column)

        node.append("g")
            .attr("transform", `translate(${marginPx.left},${marginPx.top})`)
            .attr("fill", "none")
            .attr("stroke", d3.color(color).darker(1).hex())
            .attr("stroke-linejoin", "round")
            .selectAll()
            .data(contours)
            .join("path")
            .attr("stroke-width", (d, i) => i % 5 ? 0.25 : 1)
            .attr("d", d3.geoPath(d3.geoTransform({
                // point(x, y) {
                //     this.stream.point(xAxis(x * xStd), yAxis(y * yStd))
                // }
            })));
        // for (let colI = 0; colI < columns.length; colI++) {
        //     node.append("circle").attr("cx",300).attr("cy", ).attr("r", 6).style("fill", "#69b3a2")
        //     node.append("circle").attr("cx",300).attr("cy",60).attr("r", 6).style("fill", "#404080")
        //     svg.append("text").attr("x", 320).attr("y", 30).text("variable A").style("font-size", "15px").attr("alignment-baseline","middle")
        //     svg.append("text").attr("x", 320).attr("y", 60).text("variable B").style("font-size", "15px").attr("alignment-baseline","middle")
        // }
        console.log(`KDE plot rendered in ${performance.now() - startTime}ms`)
    }

    onMount(async () => {
        await 0
        rerender()
    })

    $: {
        if (settings && data && svg) {
            Promise.resolve().then(rerender)
        }
    }
</script>

<style>

</style>

<div class="histogram-plot">
    <svg width={widthPx} height={heightPx} viewBox="0 0 {widthPx} {heightPx}">
        <g bind:this={svg}></g>
        <!-- <g>
        {#each settings?.variables ?? [] as v, i}
        {#if v.label != null}
            <circle cx={widthPx - marginPx.right - 100} cy={marginPx.top + 30 * i + 6} r="6" style="fill: {d3.color(colors[i]).brighter(0.5).hex() ?? '#000000'}; fill-opacity: 0.6; stroke: {d3.color(colors[i]).darker(0.5).hex()}"/>
            <text x={widthPx - marginPx.right - 90} y={marginPx.top + 30 * i + 12} text-anchor="start" alignment-baseline="middle">
                <title>{v.column} {#if v.filterSql}WHERE {v.filterSql}{/if}</title>
                {v.label}
            </text>
        {/if}
        {/each}
        </g> -->
    </svg>
</div>
