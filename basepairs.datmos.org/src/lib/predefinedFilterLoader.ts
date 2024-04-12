import { assetBaseUri } from "./dbInstance";
import { defaultFilter, type NucleotideFilterModel, type NumRange } from "./dbModels";
import { Lazy } from "./lazy";
import metadata from "./metadata";

export type FilterLimits = { [pairType: string]: { [column: string]: [number | null, number | null] } }

async function loadFilterCSV(url) {
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

        const family = line[columns.family]
        const bases = line[columns.bases]

        if (limits[`${family}-${bases}`] != null) {
            console.warn(`Pair type ${family}-${bases} already defined`)
            continue
        }
        limits[`${family}-${bases}`] = {}

        for (const col of limitColumns) {
            const min = line[columns[`${col}_min`]]
            const max = line[columns[`${col}_max`]]
            limits[`${family}-${bases}`][col] = [!min ? null : Number(min), !max ? null : Number(max)]
        }
    }
    return limits
}

export const defaultFilterLimits = new Lazy(async () => loadFilterCSV(new URL("filters/fr3d-data-boundaries.csv", assetBaseUri).href))

export function toNtFilter(allLimits: FilterLimits, pairType: string, baseFilter: NucleotideFilterModel | null | undefined): NucleotideFilterModel {

    const filter = baseFilter ? { ...baseFilter } : defaultFilter()

    const limits = allLimits[pairType]
    // const meta = metadata.find(m => m.pair_type.join("-").toLowerCase() == pairType.toLowerCase())
    if (!limits) {
        console.warn(`No limits for pair type ${pairType}`)
        return filter
    }

    for (const [column, [min, max]] of Object.entries(limits)) {
        const range: NumRange | undefined = min == null && max == null ? undefined : { min, max }
        // extend the range by 0.01
        if (range && range.min != null && range.max != null && range.min > range.max) {
            // inverted range
            range.min += 0.02
            range.max -= 0.02
        } else {
            if (range?.min != null)
                range.min -= 0.02
            if (range?.max != null)
                range.max += 0.02
        }
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
