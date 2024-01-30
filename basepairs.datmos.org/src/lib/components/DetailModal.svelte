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

    code {
        color: black;
        background-color: transparent;
    }

</style>
<div>
    <div class="imgpane">
        <img src={imageUrl} alt='x' />
        {#if videoUrl}
            <video src={videoUrl} autoplay muted loop controls />
        {/if}
    </div>
    <div>
        <h4>PyMol script<span style="font-size: 1rem; font-weight: 400"> &nbsp;&nbsp;- copy and paste into PyMol command line</span></h4>
        <pre>{generatePymolScript(pair.id).join("\n")}</pre>
    </div>
    <div style="display: flex; flex-direction: row;justify-content: space-evenly">
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
                    <td>{(hb.donorAngle).toFixed(0)}°</td>
                {/each}
            </tr>
            <tr>
                <th>Acceptor angle</th>
                {#each pair.hbonds as hb}
                    <td>{(hb.acceptorAngle).toFixed(0)}°</td>
                {/each}
            </tr>
        </table>
    {/if}
    {#if pair.id.nt1.pdbid}
    <div>
        <!-- <h4>From structure</h4> -->
        <div class='media' style="max-width: 600px">
            <div class="media-left">
                <a href="https://www.rcsb.org/structure/{pair.id.nt1.pdbid}">
                <code style="font-size: 3rem">{pair.id.nt1.pdbid}</code>
                </a>
            </div>
            <div class="media-content">
            <div class="content">
                <p>
                <strong>{pair.originalRow?.structure_method ?? ''}</strong> <small> at {pair.originalRow?.resolution ?? '?'} Å</small> <small>(published {pair.originalRow?.deposition_date ? new Date(pair.originalRow.deposition_date).toLocaleDateString('en', {month: 'long', day: 'numeric',  year: 'numeric'}) : ''})</small>
                <br>
                {pair.originalRow?.structure_name ?? ''}
                </p>
            </div>
            </div>
        </div>
    </div>
    {/if}
    <!-- <table class="table is-narrow is-striped" style="width: fit-content">
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
    </table> -->
    </div>
    <div>
        <table class="table is-narrow is-striped" style="width: fit-content">
            
        {#each Object.entries(pair?.originalRow ?? {}) as [k, v]}
            <tr>
                <td><b><code>{k}</code></b></td>
                <td>{typeof v == "number" ? v.toPrecision(7) : v == null ? "" : "" + v}</td>
            </tr>
        {/each}
        </table>
        <!-- <pre>{JSON.stringify(pair, null, 2)}</pre> -->
    </div>
</div>
