<script lang="ts">
	import { getColumnLabel, hideColumn, longNucleotideNames, type NucleotideFilterModel, type NumRange } from "$lib/dbModels";
	import metadata from "$lib/metadata";
	import type { NucleotideId, PairId, PairingInfo } from "$lib/pairing";
    import * as filterLoader from '$lib/predefinedFilterLoader'
	import _ from "lodash";

    export let imageUrl: string | undefined
    export let rotImageUrl: string | undefined
    export let videoUrl: string | undefined
    export let pair: PairingInfo
    export let pairType: string
    export let filter: NucleotideFilterModel | undefined = undefined
    let realFilter: NucleotideFilterModel | undefined = undefined
    $: {
        realFilter = filter
        if (filter.datasource == "allcontacts-boundaries-f") {
            filterLoader.defaultFilterLimits.value.then(v => {
                realFilter = filterLoader.addHBondLengthLimits(pairType, filterLoader.toNtFilter(v, pairType, {...filter }))
            })
        }
    }

    let videoError = false
    $: { videoUrl; videoError = false }

    function getRange(column: string) : NumRange | undefined {
        if (!realFilter) return undefined
        if (realFilter.other_column_range && column in realFilter.other_column_range) {
            return realFilter.other_column_range[column]
        }
        if (column in realFilter) {
            return realFilter[column]
        }

        if (/^C1_C1_(yaw|pitch|roll)\d*/.test(column)) {
            return realFilter[column.slice("C1_C1_".length)]
        }
    }

    function isOutOfRange(value: number | null | undefined, range: NumRange | string | null | undefined) {
        if (value == null)
            return null

        if (typeof range == "string") {
            range = getRange(range)
        }

        if (range?.min == null && range?.max == null)
            return null

        if (range.min != null && range.max != null && range.min > range.max) {
            // range in modular arithmetic (angles -180..180)
            return value > range.min || value < range.max
        }

        if (range?.min != null && value < range.min)
            return false
        if (range?.max != null && value > range.max)
            return false
        return true
    }

    function getRangeValueTitle(value: number | null | undefined, range: NumRange | string | null | undefined) {
        if (value == null || range == null)
            return ""
        if (typeof range == "string") {
            range = getRange(range)
        }

        if (range?.min == null && range?.max == null)
            return null

        const inRange = isOutOfRange(value, range)

        if (inRange === true) {
            return `Value ${value?.toFixed(2)} is within range ${range.min?.toFixed(2) ?? '-∞'}..${range.max?.toFixed(2) ?? '∞'}`
        }
        if (inRange === false) {
            return `Value ${value?.toFixed(2)} is outside range ${range.min?.toFixed(2) ?? '-∞'}..${range.max?.toFixed(2) ?? '∞'}`
        }
        return "???"
    }

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
            String(s.nt1?.chain) != String(s.nt2?.chain) ?
                `${s.nt1?.pdbid} and (chain ${s.nt1.chain} and ${resSele(s.nt1)} or chain ${s.nt2.chain} and ${resSele(s.nt2)})` :
            s.nt1?.altloc || s.nt1?.inscode || s.nt2?.altloc || s.nt2?.inscode || s.nt1.resnum < 0 || s.nt2.resnum < 0 ?
                `${s.nt1?.pdbid} and chain ${s.nt1?.chain} and (${resSele(s.nt1)} or ${resSele(s.nt2)})` :
                `${s.nt1?.pdbid} and chain ${s.nt1?.chain} and resi ${s.nt1.resnum}+${s.nt2.resnum}`
        script.push(`select pair, ${pairSelection}`)
        script.push(`show sticks, %pair`)
        script.push(`orient %pair`)
        script.push(`hide everything, not %pair`)

        return script
    }

    const columnBlacklist = [ "pdbid", "model", "chain1", "chain2", "res1", "res2", "nr1", "nr2", "alt1", "alt2", "ins1", "ins2" ]

    function getTableRows(tuple: object | null) {
        if (!tuple)  return []
        const meta = metadata.find(m => m.pair_type[0].toUpperCase() == pair.id.pairingType[0].toUpperCase() && m.pair_type[1] == pair.id.pairingType[1])

        return Object.entries(tuple).filter(([colName, _]) => !hideColumn(colName, meta) && !columnBlacklist.includes(colName)).map(([colName, value]) => {
            const [ label, tooltip ] = getColumnLabel(colName, meta) ?? [ null, null ]
            return { colName, label, tooltip, value }
        })
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

    .filter-pass {
        color: #006d09;
        font-weight: bold;
    }

    .filter-fail {
        color: rgb(172, 0, 0);
        font-weight: bold;
        text-decoration: underline;
    }

</style>
<div>
    <div class="imgpane">
        <img src={imageUrl} alt='x' />
        {#if videoUrl && !videoError}
            <video src={videoUrl} autoplay muted loop controls on:error={() => { videoError = true }} />
        {:else if rotImageUrl}
            <img src={rotImageUrl} alt='' />
        {/if}
    </div>
    <div>
        <h4>PyMol script<span style="font-size: 1rem; font-weight: 400"> &nbsp;&nbsp;- copy and paste into PyMol command line</span></h4>
        <pre>{generatePymolScript(pair.id).join("\n")}</pre>
    </div>
    <div style="display: flex; flex-direction: row; justify-content: space-evenly; gap: 2rem">
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
                {#each pair.hbonds as hb, i}
                    {@const range = realFilter?.bond_length?.[i]}
                    <td class:filter-pass={isOutOfRange(hb.length, range) === true}
                        class:filter-fail={isOutOfRange(hb.length, range) === false}
                        title={getRangeValueTitle(hb.length, range)}>
                        {hb.length == null ? 'NULL' : hb.length?.toFixed(2) + " Å"}</td>
                {/each}
            </tr>
            <tr>
                <th>Donor angle</th>
                {#each pair.hbonds as hb, i}
                    {@const range = realFilter?.bond_donor_angle?.[i]}
                    <td class:filter-pass={isOutOfRange(hb.donorAngle, range) === true}
                        class:filter-fail={isOutOfRange(hb.donorAngle, range) === false}
                        title={getRangeValueTitle(hb.donorAngle, range)}>
                        {hb.donorAngle == null ? 'NULL' : hb.donorAngle?.toFixed(0)+"°"}
                    </td>
                {/each}
            </tr>
            <tr>
                <th>Acceptor angle</th>
                {#each pair.hbonds as hb, i}
                    {@const range = realFilter?.bond_acceptor_angle?.[i]}
                    <td class:filter-pass={isOutOfRange(hb.acceptorAngle, range) === true}
                        class:filter-fail={isOutOfRange(hb.acceptorAngle, range) === false}
                        title={getRangeValueTitle(hb.acceptorAngle, range)}>
                        {hb.acceptorAngle == null ? 'NULL' : hb.acceptorAngle?.toFixed(0)+"°"}
                    </td>
                {/each}
            </tr>
            <tr>
                <th>Angle to left plane</th>
                {#each pair.hbonds as hb, i}
                    {@const range = realFilter?.bond_plane_angle1?.[i]}
                    <td class:filter-pass={isOutOfRange(hb.OOPA1, range) === true}
                        class:filter-fail={isOutOfRange(hb.OOPA1, range) === false}
                        title={getRangeValueTitle(hb.OOPA1, range)}>
                        {hb.OOPA1 == null ? 'NULL' : hb.OOPA1?.toFixed(0)+"°"}</td>
                {/each}
            </tr>
            <tr>
                <th>Angle to right plane</th>
                {#each pair.hbonds as hb, i}
                    {@const range = realFilter?.bond_plane_angle2?.[i]}
                    <td class:filter-pass={isOutOfRange(hb.OOPA2, range) === true}
                        class:filter-fail={isOutOfRange(hb.OOPA2, range) === false}
                        title={getRangeValueTitle(hb.OOPA2, range)}>
                        {hb.OOPA2 == null ? 'NULL' : hb.OOPA2?.toFixed(0)+"°"}
                    </td>
                {/each}
            </tr>
        </table>
    {/if}
    {#if pair.id.nt1.pdbid}
    <div>
        <p>
            {_.capitalize(longNucleotideNames[pair.id.nt1.resname] ?? pair.id.nt1.resname)} <strong>{pair.id.nt1.resnum}{pair.id.nt1.inscode ? '.' + pair.id.nt1.inscode : ''}</strong> in chain <strong>{pair.id.nt1.chain}</strong>
            forms {pair.id.pairingType[0]} basepair with
            {_.capitalize(longNucleotideNames[pair.id.nt2.resname] ?? pair.id.nt2.resname)} <strong>{pair.id.nt2.resnum}{pair.id.nt2.inscode ? '.' + pair.id.nt2.inscode : ''}</strong> in chain <strong>{pair.id.nt2.chain}</strong>
        </p>
        <h5 style="margin-bottom: 0px">From structure{[null, undefined, '', 0, 1, '1', '0'].includes(pair.id.nt1.model) ? '' : ` (model ${pair.id.nt1.model})`}</h5>
        <div class='media' style="max-width: 600px">
            <div class="media-left">
                <a href="https://www.rcsb.org/structure/{pair.id.nt1.pdbid}">
                <code style="font-size: 3rem">{pair.id.nt1.pdbid}</code>
                </a>
            </div>
            <div class="media-content">
            <div class="content">
                <p>
                <strong>{pair.originalRow?.structure_method ?? ''}</strong> <small> at {pair.originalRow?.resolution?.toFixed(2) ?? '?'} Å</small> <small>(published {pair.originalRow?.deposition_date ? new Date(pair.originalRow.deposition_date).toLocaleDateString('en', {month: 'long', day: 'numeric',  year: 'numeric'}) : ''})</small>
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
            
        <tbody>
        {#each getTableRows(pair?.originalRow) as r}
            {@const filterRange = getRange(r.colName)}
            {@const val = r.value}
            <tr>
                <td><b><code>{r.colName}</code></b></td>
                <td>{r.label ?? ''}</td>
                <td colspan={r.colName == 'structure_name' ? 2 : 1}
                    style="font-weigth: 700; text-align: right;"
                    class:filter-pass={isOutOfRange(val, filterRange) === true}
                    class:filter-fail={(val == null && filterRange != null) ||  isOutOfRange(val, filterRange) === false}
                    title={getRangeValueTitle(val, filterRange)}
                    data-type={typeof val}>
                    {typeof val == "number" ? val.toFixed(3) : val == null ? (filterRange != null ? "NULL" : "") : "" + val}
                </td>
                <td><i>{r.tooltip ?? ''}</i></td>
            </tr>
        {/each}
        </tbody>
        </table>
        <!-- <pre>{JSON.stringify(pair, null, 2)}</pre> -->
    </div>
</div>
