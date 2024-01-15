<script lang="ts">
	import type { NucleotideId, PairId, PairingInfo } from "$lib/pairing";

    export let imageUrl: string
    export let videoUrl: string
    export let pair: PairingInfo


    function resSele(r: NucleotideId) {
        let nt = String(r.resnum).replace("-", "\\-")
        if (r.inscode) {
            nt += String(r.inscode)
        }

        const alt = r.altloc ? ` alt ${r.altloc}` : ""
        return `resi ${nt}${alt}`
    }
    
    function generatePymolScript(s: PairId): string[] {
        const script = []
        script.push(`fetch ${s.nt1?.pdbid}`)
        const pairSelection =
            String(s.nt1?.chain) == String(s.nt2?.chain) ?
                `${s.nt1?.pdbid} and chain ${s.nt1?.chain} and (${resSele(s.nt1)} or ${resSele(s.nt2)})` :
                `${s.nt1?.pdbid} and (chain ${s.nt1.chain} and ${resSele(s.nt1)} or chain ${s.nt2.chain} and ${resSele(s.nt2)})`;
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
    table th {
        white-space: nowrap;
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
    <div style="display: flex; flex-direction: row">
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
    <table class="table is-narrow is-striped" style="width: fit-content">
        <tr>
            <th>Structure ID</th>
            <td>{pair.id.nt1.pdbid}</td>
        </tr>
        <tr>
            <th>Structure Name</th>
            <td>{pair.originalRow?.structure_name ?? ''}</td>
        </tr>
        <tr>
            <th>Structure Method</th>
            <td>{pair.originalRow?.structure_method ?? ''}</td>
        </tr>
        <tr>
            <th>Resolution </th>
            <td>{pair.originalRow?.resolution ?? '?'} Å</td>
        </tr>
        <tr>
            <th>Deposition date</th>
            <td>{pair.originalRow?.deposition_date ? new Date(pair.originalRow.deposition_date).toLocaleDateString('en', {month: 'long', day: 'numeric',  year: 'numeric'}) : ''}</td>
        </tr>
    </table>
    </div>
    <div>
        {#each Object.entries(pair?.originalRow ?? {}) as [k, v]}
            <div><b>{k}</b>: {v}</div>
        {/each}
        <!-- <pre>{JSON.stringify(pair, null, 2)}</pre> -->
    </div>
</div>
