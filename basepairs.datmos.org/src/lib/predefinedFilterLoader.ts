import config from "./config";
import { assetBaseUri } from "./dbInstance";
import { defaultFilter, type NucleotideFilterModel, type NumRange } from "./dbModels";
import { Lazy } from "./lazy";
import metadata from "./metadata";

export type FilterLimits = { [pairType: string]: { [column: string]: [number | null, number | null] } }

async function loadFilterCSV(url) {
    url = new URL(url, assetBaseUri).href
    const response = await fetch(url)
    const lines = (await response.text()).split("\n").map(line => line.trim())

    const columns: {[col: string]: number} = Object.fromEntries(lines[0].split(",").map((column, index) => [column, index]))

    const limitColumns = [...new Set(Object.keys(columns).filter(column => column.endsWith("_min") || column.endsWith("_max")).map(column => column.slice(0, -4)))]

    const limits: FilterLimits = {}

    for (let i = 1; i < lines.length; i++) {
        const line = lines[i].split(",").map(value => value.trim())
        if (line.length == 1) {
            continue
        }

        const key = `${line[columns.family]}-${line[columns.bases]}`.toLowerCase()

        if (limits[key] != null) {
            console.warn(`Pair type key already defined`)
            continue
        }
        limits[key] = {}

        for (const col of limitColumns) {
            const min = line[columns[`${col}_min`]]
            const max = line[columns[`${col}_max`]]
            limits[key][col] = [!min ? null : Number(min), !max ? null : Number(max)]
        }
    }
    return limits
}

export const filterLimits = new Lazy(async () => {
    const e = Object.entries(config.parameterBoundariesUrls)
    const f = await Promise.all(e.map(x => loadFilterCSV(x[1]).catch(err => {
        console.error(`Failed to load filter limits ${x[0]} from ${x[1]}`, err)
        return null
    })))
    return Object.fromEntries(e.map((x, i) => [x[0], f[i]]))
})

export const defaultFilterLimits = new Lazy(() => filterLimits.value.then(l => {
    const result = config.defaultBoundaries.find(b => l[b] != null)
    if (!result) {
        throw new Error(`No default filter limits could be loaded`)
    }
    return { id: result, limits: l[result] }
}))

export function addHBondLengthLimits(pairType: string, extendBy: number, baseFilter: NucleotideFilterModel) {
    pairType = pairType.toLowerCase()
    const m = metadata.find(m => m.pair_type.join("-").toLowerCase() == pairType)
    const hbLengthLimits = m.atoms.map(([_, a, b, __], i) =>
        (a.includes("C") || b.includes("C") || b.includes("O2'") || a.includes("O2'") || i >= 2 ? 4.2 : 4) + extendBy)
    return {
        ...baseFilter,
        min_bond_length: baseFilter.min_bond_length ?? {max: 3.8},
        bond_length: hbLengthLimits.map((l, i) => baseFilter.bond_length[i] || ({ min: null, max: l }))
    }
}

export function toNtFilter(allLimits: FilterLimits, extendBy: number, pairType: string, baseFilter: NucleotideFilterModel | null | undefined): NucleotideFilterModel {

    const filter = baseFilter ? structuredClone(baseFilter) : defaultFilter()

    const limits = allLimits[pairType.toLowerCase()]
    // const meta = metadata.find(m => m.pair_type.join("-").toLowerCase() == pairType.toLowerCase())
    if (!limits) {
        console.warn(`No limits for pair type ${pairType}`)
        return filter
    }

    for (const [column, [min, max]] of Object.entries(limits)) {
        const range: NumRange | undefined = min == null && max == null ? undefined : { min, max }
        // extend the range by 0.01
        if (range?.min != null)
            range.min -= extendBy
        if (range?.max != null)
            range.max += extendBy
        if (/^hb_\d/.test(column)) {
            const [_, hbIndex, hbParam] = /hb_(\d+)_(.+)/.exec(column)
            const hb = Number(hbIndex)

            const arr = hbParam == "length" ? filter.bond_length :
                        hbParam == "AA" || hbParam == "acceptor_angle" ? filter.bond_acceptor_angle :
                        hbParam == "DA" || hbParam == "donor_angle" ? filter.bond_donor_angle :
                        hbParam == "OOPA1" || hbParam == "plane_angle1" ? filter.bond_plane_angle1 :
                        hbParam == "OOPA2" || hbParam == "plane_angle2" ? filter.bond_plane_angle2 :
                        null

            if (arr != null) {
                while (arr.length <= hb) {
                    arr.push(undefined)
                }
                arr[hb] = range
                continue
            }
        }

        if (column == "min_bond_length") {
            filter.min_bond_length = range
            continue
        }

        if (column == "coplanarity_angle") {
            filter.coplanarity_angle = range
            continue
        }
        if (column == "coplanarity_edge_angle1") {
            filter.coplanarity_edge_angle1 = range
            continue
        }
        if (column == "coplanarity_edge_angle2") {
            filter.coplanarity_edge_angle2 = range
            continue
        }
        if (column == "coplanarity_shift1") {
            filter.coplanarity_shift1 = range
            continue
        }
        if (column == "coplanarity_shift2") {
            filter.coplanarity_shift2 = range
            continue
        }
        if (column == "yaw1" || column == "C1_C1_yaw1") {
            filter.yaw1 = range
            continue
        }
        if (column == "pitch1" || column == "C1_C1_pitch1") {
            filter.pitch1 = range
            continue
        }
        if (column == "roll1" || column == "C1_C1_roll1") {
            filter.roll1 = range
            continue
        }
        if (column == "yaw2" || column == "C1_C1_yaw2") {
            filter.yaw2 = range
            continue
        }
        if (column == "pitch2" || column == "C1_C1_pitch2") {
            filter.pitch2 = range
            continue
        }
        if (column == "roll2" || column == "C1_C1_roll2") {
            filter.roll2 = range
            continue
        }

        filter.other_column_range[column] = range
    }

    return filter

}
