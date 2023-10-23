<script lang="ts">
	import type { NucleotideId, PairId, PairingInfo } from "$lib/pairing";

    export let imageUrl: string
    export let videoUrl: string
    export let pair: PairingInfo


    function residueSelection(r: NucleotideId) {
        const chain = String(r.chain)
        let nt = String(r.resnum).replace("-", "\\-")
        if (r.inscode) {
            nt += String(r.inscode)
        }

        const alt = r.altloc ? ` alt ${r.altloc}` : ""
        return `(chain ${chain} and resi ${nt}${alt})`
    }
    
    function generatePymolScript(s: PairId): string[] {
        const script = []
        script.push(`fetch ${s.nt1?.pdbid}`)
        const pairSelection = `${s.nt1?.pdbid} and (${residueSelection(s.nt1)} or ${residueSelection(s.nt2)})`
        script.push(`select pair, ${pairSelection}`)
        script.push(`show sticks, %pair`)
        script.push(`orient %pair`)
        script.push(`hide everything, not %pair`)

        return script
    }
</script>

<style>
    .imgpane {
        display: flex;
        align-items: center;
    }
    .imgpane > * {
        flex-grow: 1;
        width: 40vw;
    }


</style>
<div>
    <div class="imgpane">
        <img src={imageUrl} alt='x' />
        {#if videoUrl}
            <video src={videoUrl} autoplay muted loop />
        {/if}
    </div>
    <div>
        <h4>PyMol script<span style="font-size: 1rem; font-weight: 400"> &nbsp;&nbsp;- copy and paste into PyMol command line</span></h4>
        <pre>{generatePymolScript(pair.id).join("\n")}</pre>
    </div>
    <div>
    {#if pair.hbonds}
        <table class="table is-narrow is-striped" style="width: fit-content">
            <tr>
                <th></th>
                {#each pair.hbonds as hb, ix}
                    <th>{hb.label ?? `H-Bond ${ix}`}</th>
                {/each}
            </tr>
            <tr>
                <th>Length</th>
                {#each pair.hbonds as hb}
                    <td>{hb.length.toFixed(2)} Å</td>
                {/each}
            </tr>
            <tr>
                <th>Donor angle</th>
                {#each pair.hbonds as hb}
                    <td>{(hb.donorAngle * 180 / Math.PI).toFixed(0)}°</td>
                {/each}
            </tr>
            <tr>
                <th>Acceptor angle</th>
                {#each pair.hbonds as hb}
                    <td>{(hb.acceptorAngle * 180 / Math.PI).toFixed(0)}°</td>
                {/each}
            </tr>
        </table>
    {/if}
    </div>
    <div>
        <!-- <pre>{JSON.stringify(pair, null, 2)}</pre> -->
    </div>
</div>
