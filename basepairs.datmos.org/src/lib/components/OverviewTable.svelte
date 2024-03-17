A<script lang="ts">
    import { pairFamilies, pairTypes } from "$lib/dbInstance";
    import metadata from "$lib/metadata";
	import Pairimage from "./pairimage.svelte";
    import { imgDir } from "$lib/dbInstance";
    const allFamilies = [ "cWW", "tWW", "cWH", "tWH", "cWS", "tWS", "cHH", "tHH", "cHS", "tHS", "cSS", "tSS" ];
    export let families: string[] | undefined = undefined

    const bases = [ "A", "G", "C", "U" ];

    function familyLabel(family: string) {
        const n = family.startsWith("n")
        if (n) family = family.slice(1);
        const cis = family.startsWith("c");
        const mapping = { "W": "Watson-Crick", "H": "Hoogsteen", "S": "Sugar", "B": "Bifurcated" };
        const ix = allFamilies.indexOf(family);
        return `${ix+1}: ${n ? "nearly " : ""} ${cis ? "cis" : "trans"} ${mapping[family[1]]} / ${mapping[family[2]]}`;
    }
    // function getUrl(meta: ) {
    //     if (!pair?.nt1.pdbid || !pair?.nt2.pdbid || pair.nt1.chain == null || pair.nt2.chain == null || pair.nt1.resnum == null || pair.nt2.resnum == null)
    //         return undefined
    //     const imgName = `${pair.nt1.chain}_${pair.nt1.resnum}${pair.nt1.inscode??''}${pair.nt1.altloc??''}-${pair.nt2.chain}_${pair.nt2.resnum}${pair.nt2.inscode??''}${pair.nt2.altloc??''}`
    //     return `${imgDir}/${pair.nt1.pdbid}/${imgName}${attachement}`
    // }
</script>

<style>
    .LW-table {
        display: flex;
        flex: 0 0 auto;
        flex-direction: row;
        justify-content: space-between;
        flex-wrap: wrap;
        align-items: center;
        margin: 1rem 0;
        padding-left: 3rem;
        width: 100vw;
        --cell-height: calc(var(--cell-width) / 16 * 9 + 2rem);
    }
    @media (max-width: 1000px) {
        .LW-table {
            --cell-width: calc((100vw - 6rem) / 2);
        }
    }
    @media (min-width: 1000px) {
        .LW-table {
            --cell-width: calc((100vw - 6rem) / 4);
        }
    }

    .LW-table-column {
        display: flex;
        flex-direction: column;
        justify-content: space-between;
        align-items: center;
    }
    .LW-table-column:not(.left-labels) {
        width: var(--cell-width);
    }
    .LW-table-column.left-labels {
        margin-left: -2rem;
        width: 2rem;
    }
    .LW-table-column.left-labels > div:not(.LW-table-col-head) {
        padding-top: calc((var(--cell-height) - 2rem) / 2);
        font-size: 2rem;
    }


    .LW-table-col-head {
        height: 3rem;
        font-size: 2rem;
    }
    .LW-table-column > div:not(.LW-table-col-head) {
        height: var(--cell-height);
        text-align: center;
        display: flex;
        flex-direction: column;
        justify-content: center;
        align-items: center;
    }

    h2 {
        font-size: 2rem;
        text-align: center;
        margin-top: 0rem;
        margin-bottom: 0rem;
        border-top: solid 1px black;
    }
</style>

<div>
    {#each (families ?? allFamilies) as family}
    <h2>{familyLabel(family)}</h2>
    <div data-LW-family={family} class="LW-table">
        <div class="LW-table-column left-labels">
            <div class="LW-table-col-head"></div>
            {#each bases as base}
                <div>{base}</div>
            {/each}
        </div>
        {#each bases as base1}
            <div class="LW-table-column">
                <div class="LW-table-col-head">{base1}</div>
                {#each bases as base2}
                    {@const meta = metadata.find(m => m.pair_type[0].toUpperCase() == family.toUpperCase() && m.pair_type[1] == `${base1}-${base2}`)}
                    {@const biggestStat = meta?.statistics.reduce((prev, x) => (prev && prev.count > x.count ? prev : x), null)}
                    {@const url = biggestStat?.nicest_bp == null ? null : `${imgDir}/${biggestStat.nicest_bp[0]}/${biggestStat.nicest_bp[2]}_${biggestStat.nicest_bp[4]}${biggestStat.nicest_bp[5]??''}${biggestStat.nicest_bp[6]??''}-${biggestStat.nicest_bp[7]}_${biggestStat.nicest_bp[9]}${biggestStat.nicest_bp[10]??''}${biggestStat.nicest_bp[11]??''}.png`}
                    <div class="LW-table-cell" style="text-align: center">
                        {#if url != null}
                        <!-- onClick={() => location.hash = `#${family}-${base1}-${base2}`} -->
                        <Pairimage
                            parentSize={true}
                            url={url}
                            videoUrl={undefined}
                            allowHoverVideo={false}
                            linkText="show statistics + exemplars"
                            linkUrl={`#${family}-${base1}-${base2}`}
                            onClick={() => true}
                            pair={{
                                id: {
                                    nt1: {pdbid: biggestStat.nicest_bp[0], model: biggestStat.nicest_bp[1], chain: biggestStat.nicest_bp[2], resname: biggestStat.nicest_bp[3], resnum: biggestStat.nicest_bp[4], altloc: biggestStat.nicest_bp[5], inscode: biggestStat.nicest_bp[6]},
                                    nt2: {pdbid: biggestStat.nicest_bp[0], model: biggestStat.nicest_bp[1], chain: biggestStat.nicest_bp[7], resname: biggestStat.nicest_bp[8], resnum: biggestStat.nicest_bp[9], altloc: biggestStat.nicest_bp[10], inscode: biggestStat.nicest_bp[11]},
                                    pairingType: [ family, `${base1}-${base2}` ]
                                }
                            }} />
                        {:else if
                            ['ww', 'hh'].includes(family.toLowerCase().replace(/^n?[ct]/, '')) &&
                            metadata.find(m => m.pair_type[0].toUpperCase() == family.toUpperCase() && m.pair_type[1] == `${base2}-${base1}`)}
                            <div style="font-size: 1.5rem; color: rgb(170, 170, 170)">same as <b>{base2}-{base1}</b></div>
                        {:else}
                            <div style="font-size: 1.5rem; color: rgb(255, 170, 170)">not defined</div>
                        {/if}
                    </div>
                {/each}
            </div>
        {/each}
    </div>
    {/each}
</div>
