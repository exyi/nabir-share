import type metadataModule from './metadata'
import _ from 'lodash'
import type { PairingInfo } from './pairing'
import config from '$lib/config'


export type NumRange = {
    min?: number
    max?: number
}
export type NucleotideFilterModel = {
    sql?: string
    bond_length: (NumRange)[]
    bond_donor_angle: (NumRange)[]
    bond_acceptor_angle: (NumRange)[]
    bond_plane_angle1: (NumRange)[]
    bond_plane_angle2: (NumRange)[]
    other_column_range?: {
        [column: string]: NumRange
    }
    sql_conditions?: string[]
    coplanarity_angle?: NumRange
    yaw1?: NumRange
    yaw2?: NumRange
    pitch1?: NumRange
    pitch2?: NumRange
    roll1?: NumRange
    roll2?: NumRange
    coplanarity_shift1?: NumRange
    coplanarity_shift2?: NumRange
    coplanarity_edge_angle1?: NumRange
    coplanarity_edge_angle2?: NumRange
    resolution?: NumRange
    dna?: true | false | undefined
    orderBy: string
    filtered: boolean
    datasource?: "fr3d" | "fr3d-f" | "fr3d-n" | "fr3d-nf" | "allcontacts-f" | "allcontacts-boundaries-f"
    includeNears?: boolean
    rotX?: boolean
}

export type StatisticsSettingsModel = {
    enabled: boolean
    panels: (StatPanelSettingsModel)[]
}
export type StatPanelSettingsModel = HistogramSettingsModel | KDE2DSettingsModel
export type HistogramSettingsModel = {
    type: "histogram"
    title?: string
    bins?: number
    variables: VariableModel[]
}
export type KDE2DSettingsModel = {
    type: "kde2d"
    title?: string
    variables: VariableModel[]
}
export type VariableModel = {
    column: string
    label: string
    tooltip?: string
    filterSql?: string
    filterId?: string
}

export type DetailModalViewModel = {
    pair: PairingInfo
    imgUrl?: string
    rotImgUrl?: string
    videoUrl?: string

}

export function defaultFilter(datasource: NucleotideFilterModel["datasource"] = config.defaultDataSource): NucleotideFilterModel {
    return { datasource, bond_acceptor_angle: [], bond_donor_angle: [], bond_length: [], bond_plane_angle1: [], bond_plane_angle2: [], filtered: true, orderBy: config.defaultOrderBy }
}

function rangeToCondition(col: string, range: NumRange | undefined | null): string[] {
    if (range == null) return []

    function almostEqual(a: number, b: number) {
        return Math.abs(a - b) < (0.001 * Math.max(Math.abs(a), Math.abs(b)))
    }

    if (range.min != null && range.max != null) {
        if (range.min > range.max) {
            // modular range
            return [`${col} NOT BETWEEN ${range.max} AND ${range.min}`]
        }
        if (range.min == range.max) {
            return [`${col} = ${range.min}`]
        }
        if (almostEqual(range.min, -range.max)) {
            return [`abs(${col}) <= ${Math.max(Math.abs(range.min), Math.abs(range.max))}`]
        }
        return [`${col} BETWEEN ${range.min} AND ${range.max}`]
    }
    if (range.min != null) {
        return [`${col} >= ${range.min}`]
    }
    if (range.max != null) {
        return [`${col} <= ${range.max}`]
    }
    return []
}

export function filterToSqlCondition(filter: NucleotideFilterModel): string[] {
    const conditions: string[] = []
    if (filter.filtered) {
        conditions.push(`jirka_approves`)
    }
    for (let i = 0; i < 3; i++) {
        conditions.push(...rangeToCondition(`hb_${i}_length`, filter.bond_length[i]))
        conditions.push(...rangeToCondition(`hb_${i}_donor_angle`, filter.bond_donor_angle[i]))
        conditions.push(...rangeToCondition(`hb_${i}_acceptor_angle`, filter.bond_acceptor_angle[i]))
        conditions.push(...rangeToCondition(`hb_${i}_OOPA1`, filter.bond_plane_angle1[i]))
        conditions.push(...rangeToCondition(`hb_${i}_OOPA2`, filter.bond_plane_angle2[i]))
    }
    conditions.push(...rangeToCondition(`resolution`, filter.resolution))
    if (filter.coplanarity_angle) {
        conditions.push(...rangeToCondition(`coplanarity_angle`, filter.coplanarity_angle))
    }
    if (filter.coplanarity_shift1) {
        conditions.push(...rangeToCondition(`coplanarity_shift1`, filter.coplanarity_shift1))
    }
    if (filter.coplanarity_shift2) {
        conditions.push(...rangeToCondition(`coplanarity_shift2`, filter.coplanarity_shift2))
    }
    if (filter.coplanarity_edge_angle1) {
        conditions.push(...rangeToCondition(`coplanarity_edge_angle1`, filter.coplanarity_edge_angle1))
    }
    if (filter.coplanarity_edge_angle2) {
        conditions.push(...rangeToCondition(`coplanarity_edge_angle2`, filter.coplanarity_edge_angle2))
    }
    if (filter.yaw1) {
        conditions.push(...rangeToCondition(`C1_C1_yaw1`, filter.yaw1))
    }
    if (filter.yaw2) {
        conditions.push(...rangeToCondition(`C1_C1_yaw2`, filter.yaw2))
    }
    if (filter.pitch1) {
        conditions.push(...rangeToCondition(`C1_C1_pitch1`, filter.pitch1))
    }
    if (filter.pitch2) {
        conditions.push(...rangeToCondition(`C1_C1_pitch2`, filter.pitch2))
    }
    if (filter.roll1) {
        conditions.push(...rangeToCondition(`C1_C1_roll1`, filter.roll1))
    }
    if (filter.roll2) {
        conditions.push(...rangeToCondition(`C1_C1_roll2`, filter.roll2))
    }

    if (filter.dna != null) {
        if (filter.dna)
            conditions.push(`(res1 LIKE 'D%' OR res2 LIKE 'D%')`)
        else
            conditions.push(`(res1 NOT LIKE 'D%' OR res2 NOT LIKE 'D%')`)
    }
    return conditions
}

export function getDataSourceTable(filter: NucleotideFilterModel) {
    const ds = filter.datasource ?? config.defaultDataSource
    if (ds == "fr3d-f") {
        return "selectedpair_f"
    } else if (ds == "fr3d") {
        return "selectedpair"
    } else if (ds == "fr3d-n") {
        return "(select * FROM selectedpair UNION ALL BY NAME SELECT * from selectedpair_n)"
    } else if (ds == "fr3d-nf") {
        return "(select * FROM selectedpair_f UNION ALL BY NAME SELECT * from selectedpair_n where jirka_approves)"
    } else if (ds == "allcontacts-f") {
        return "selectedpair_allcontacts_f"
    } else if (ds == "allcontacts-boundaries-f") {
        return "selectedpair_allcontacts_boundaries_f"
    } else {
        console.error("Unknown datasource", ds)
        return "selectedpair_allcontacts"
    }
}


function joinConditions(conditions: string[], keyword = 'AND') {
    if (conditions.length == 0) {
        return keyword.toUpperCase() != 'or' ? 'true' : 'false'
    }
    return conditions.map(c => /\b(select|or)\b/.test(c) ? `(${c})` : c).join(`\n  ${keyword} `)
}

function buildSelect(opt: {
    cols: string,
    from: string,
    where?: string[]
    limit?: number | string
    orderBy?: string
}) {
    let query = `SELECT ${opt.cols} FROM ${opt.from}`
    if (opt.where && opt.where?.length > 0) {
        query += `\nWHERE ${joinConditions(opt.where)}`
    }
    if (opt.orderBy) {
        query += `\nORDER BY ${orderToExpr(opt.orderBy)}`
    }
    if (opt.limit) {
        query += `\nLIMIT ${opt.limit}`
    }
    return query
}

export function makeSqlQuery(filter: NucleotideFilterModel, from: string, limit?: number, cols = "*") {
    return buildSelect({
        cols,
        from,
        where: filterToSqlCondition(filter),
        orderBy: filter.orderBy,
        limit
    })
}

export const orderByOptions = [
    { id: "pdbid", expr: "", label: "pdbid" },
    { id: "pdbidD", expr: "pdbid DESC, model DESC, chain1 DESC, nr1 DESC", title: "", label: "pdbid descending" },
    { id: "resolutionA", expr: "resolution NULLS LAST, pdbid, model, chain1, nr1", title: "Reported resolution of the source PDB structure", value: "resolution" },
    { id: "modedevA", expr: "mode_deviations", title: "ASCENDING - best to worst - Number of standard deviations between the H-bond parameters and the modes (peaks) calculated from Kernel Density Estimate. Use to list &quot;nicest&quot; pairs and avoid secondary modes.", label: "Deviation from KDE mode ↓" },
    { id: "modedevD", expr: "-mode_deviations", title: "DESCENDING - worst to best - Number of standard deviations between the H-bond parameters and the modes (peaks) calculated from Kernel Density Estimate. Use to list &quot;nicest&quot; pairs and avoid secondary modes.", label: "Deviation from KDE mode ↑" },
    { id: "LLD", expr: "-log_likelihood", title: "↑ DESCENDING - best to worst - Multiplied likelihoods of all H-bond parameters in their Kernel Density Estimate distribution. Use to list &quot;nicest&quot; pairs without disqualifying secondary modes.", label: "KDE likelihood ↑" },
    { id: "LLA", expr: "log_likelihood", title: "↓ ASCENDING - best to worst - Multiplied likelihoods of all H-bond parameters in their Kernel Density Estimate distribution. Use to list &quot;nicest&quot; pairs without disqualifying secondary modes.", label: "KDE likelihood ↓" },
    { id: "rmsdA", expr: "(rmsd_edge1+rmsd_edge2)", title: "↓ ASCENDING ('best first') - edge RMSD to the 'nicest' basepair", label: "Smallest Edge RMSD" },
    { id: "rmsdD", expr: "(rmsd_edge1+rmsd_edge2) DESC", title: "↑ DESCENDING ('worst first') - edge RMSD to the 'nicest' basepair", label: "Largest Edge RMSD" }
]

function orderToExpr(opt: string | null, prefix = '') {
    if (!opt) {
        return null
    }

    return prefix + (orderByOptions.find(o => o.id == opt)?.expr ?? opt)
}

export type ComparisonMode = "union" | "difference" | "missing" | "new"

export function makeDifferentialSqlQuerySingleTable(filter: NucleotideFilterModel, filterBaseline: NucleotideFilterModel, from: string, limit: number | null, mode: ComparisonMode, columns = "*") {
    const conditions = filterToSqlCondition(filter)
    const baselineConditions = filterToSqlCondition(filterBaseline)

    const commonConditions = conditions.filter(c => baselineConditions.includes(c))

    if (mode == "missing" || mode == "new") {
        // missing = is in baseline but not in current
        //         = SELECT * FROM dataset WHERE baseline_filter AND NOT (current_filter)
        const missing = mode == "new"
        const not = (mode == "missing" ? conditions : baselineConditions).filter(c => !commonConditions.includes(c))
        return buildSelect({
            cols: `${columns},
                ${missing} AS comparison_in_baseline,
                ${!missing} AS comparison_in_current`,
            from,
            where: [
                ...(mode == "missing" ? baselineConditions : conditions),
                `NOT (${joinConditions(not, 'AND')})`
            ],
            orderBy: filter.orderBy,
            limit
        })
    }


    const queryBase = buildSelect({
        cols: "*",
        from,
        where: commonConditions,
        orderBy: filter.orderBy
    })

    let query = `
        SELECT ${columns},
            coalesce(${joinConditions(conditions.filter(c => !commonConditions.includes(c)))}, FALSE) AS comparison_in_current,
            coalesce(${joinConditions(baselineConditions.filter(c => !commonConditions.includes(c)))}, FALSE) AS comparison_in_baseline
        FROM (${queryBase})
        WHERE ${mode == "difference" ? 'comparison_in_current <> comparison_in_baseline' : 'comparison_in_current OR comparison_in_baseline'}
    `
    if (limit) {
        query += `\nLIMIT ${limit}`
    }
    return query
}

export function makeDifferentialSqlQuery(queryCurrent: string, queryBaseline: string, limit?: number, order?: string, mode: ComparisonMode = "difference") {
    let query
    if (mode == "new" || mode == "missing") {
        const missing = mode == "missing"
        query = `
            SELECT *,
                ${missing} AS comparison_in_baseline,
                ${!missing} AS comparison_in_current
            FROM (${missing ? queryBaseline : queryCurrent})
            WHERE pairid NOT IN (SELECT pairid FROM (${missing ? queryCurrent : queryBaseline}))
        `
    } else {
        query = `
            SELECT DISTINCT ON (pairid)
                * EXCLUDE (comparison_in_baseline, comparison_in_current),
                bool_or(comparison_in_baseline) OVER (PARTITION BY pairid) as comparison_in_baseline,
                bool_or(comparison_in_current) OVER (PARTITION BY pairid) as comparison_in_current
            FROM (
                SELECT *, TRUE AS comparison_in_current, FALSE AS comparison_in_baseline FROM (${queryCurrent})
                UNION ALL BY NAME
                SELECT *, FALSE AS comparison_in_current, TRUE AS comparison_in_baseline FROM (${queryBaseline})
                ORDER BY comparison_in_current DESC
            )
        `
        if (mode == "difference") {
            query = `SELECT * FROM (${query})\nWHERE comparison_in_baseline <> comparison_in_current`
        } else {
            query = `SELECT * FROM (${query})` // must wrap to run ORDER BY after the SELECT DISTINCT ON (pairid)
        }
    }
    if (limit) {
        query += `\nLIMIT ${limit}`
    }
    if (order) {
        query += `\nORDER BY ${orderToExpr(order)}`
    }
    return query
}

export function aggregateTypesQuery(query: string, family = "family", res1 = "res1", res2 = "res2") {
    return `
        SELECT concat(${family}, '-', ltrim(${res1}, 'D'), '-', ltrim(${res2}, 'D')) as type,
               COUNT(*) AS count
        FROM (${query})
        GROUP BY ${family}, ltrim(${res1}, 'D'), ltrim(${res2}, 'D')
        ORDER BY COUNT(*) DESC`
}

export function aggregatePdbCountQuery(query: string, pdbid = "pdbid") {
    return `
        SELECT ${pdbid} as pdbid, COUNT(*) AS count
        FROM (${query})
        GROUP BY ${pdbid}
        ORDER BY COUNT(*) DESC`
}

export function aggregateComparisoonTypeQuery(query: string) {
    return `
        SELECT comparison_in_baseline, comparison_in_current, COUNT(*) AS count
        FROM (${query})
        GROUP BY comparison_in_baseline, comparison_in_current
        ORDER BY comparison_in_baseline, comparison_in_current`
}

export function aggregateCountsAcrossGroups(query: string, groupingSets: string[][], agg = { count: 'COUNT(*)' }) {
    const columns = []
    const colHashSet = new Set<string>()
    for (const set of groupingSets) {
        for (const col of set) {
            if (colHashSet.has(col)) continue
            columns.push(col)
            colHashSet.add(col)
        }
    }
    return `
        SELECT ${Object.entries(agg).map(([name, expr]) => `${expr} AS ${name}`).join(', ')},
                ${columns.map(c => `${c}`).join(', ')}
        FROM (${query})
        GROUP BY GROUPING SETS (${groupingSets.map(set => `(${set.join(', ')})`).join(', ')})`
}

export function aggregateBondParameters(query: string, bondColumns: string[]) {
    // const columns = bondColumns.map(c => {
    //     const m = /hb_(\d)_(\w+)/.exec(c)
    //     if (!m) return null
    //     const [_, i, name] = m
    //     return [ Number(i), name ]
    // }).filter(c => c != null)
    const columns = bondColumns.filter(c => /^hb_\d_\w+/.test(c))
    return `
        SELECT ${columns.map((c) => `
            avg(${c}) as ${c}_avg,
            min(${c}) as ${c}_min,
            max(${c}) as ${c}_max,
            median(${c}) as ${c}_median,
            stddev(${c}) as ${c}_stddev,
            count(${c}) as ${c}_nncount
            `).join(', ')},
        FROM (${query}) data
    `;

//     approx_quantile(${c}, 0.10) as ${c}_p10,
//     approx_quantile(${c}, 0.25) as ${c}_p25,
//     approx_quantile(${c}, 0.75) as ${c}_p75,
//     approx_quantile(${c}, 0.90) as ${c}_p90,
}


export function filterToUrl(filter: NucleotideFilterModel, filterBaseline: NucleotideFilterModel | null | undefined = undefined, mode = 'ranges') {
    const params = new URLSearchParams()
    addFilterParams(params, filter, mode, '')
    if (filterBaseline != null) {
        addFilterParams(params, filterBaseline, filterBaseline.sql ? "sql" : "ranges", 'baseline_')
    }
    return params
}

export function hasFilters(filter: NucleotideFilterModel, mode: string) {
    const params = new URLSearchParams()
    addFilterParams(params, filter, mode)
    params.delete('ds')
    params.delete('order')
    return params.size > 0
}

export function addFilterParams(params: URLSearchParams, filter: NucleotideFilterModel, mode: string, prefix='') {
    if (mode == 'sql') {
        add('sql', filter.sql)
        return;
    }

    if (prefix != '' || filter.datasource != null && filter.datasource != config.defaultDataSource || filter.dna != null) {
        add('ds', (filter.datasource ?? config.defaultDataSource) + (filter.dna == true ? 'D' : filter.dna == false ? 'R' : ''))
    }
    for (let i = 0; i < 3; i++) {
        addMaybe(`hb${i}_L`, range(filter.bond_length[i]))
        addMaybe(`hb${i}_DA`, range(filter.bond_donor_angle[i]))
        addMaybe(`hb${i}_AA`, range(filter.bond_acceptor_angle[i]))
        addMaybe(`hb${i}_OOPA1`, range(filter.bond_plane_angle1[i]))
        addMaybe(`hb${i}_OOPA2`, range(filter.bond_plane_angle2[i]))
    }
    addMaybe(`coplanarity_a`, range(filter.coplanarity_angle))
    addMaybe(`coplanarity_edge_angle1`, range(filter.coplanarity_edge_angle1))
    addMaybe(`coplanarity_edge_angle2`, range(filter.coplanarity_edge_angle2))
    addMaybe(`coplanarity_shift1`, range(filter.coplanarity_shift1))
    addMaybe(`coplanarity_shift2`, range(filter.coplanarity_shift2))
    addMaybe(`yaw1`, range(filter.yaw1))
    addMaybe(`yaw2`, range(filter.yaw2))
    addMaybe(`pitch1`, range(filter.pitch1))
    addMaybe(`pitch2`, range(filter.pitch2))
    addMaybe(`roll1`, range(filter.roll1))
    addMaybe(`roll2`, range(filter.roll2))
    for (const c of filter.sql_conditions ?? []) {
        add('condition', c)
    }
    for (const [column, r] of Object.entries(filter.other_column_range ?? {})) {
        addMaybe("range_" + column, range(r))
    }
    addMaybe(`resol`, range(filter.resolution))
    if (filter.orderBy != null && filter.orderBy != config.defaultOrderBy) {
        addMaybe(`order`, filter.orderBy)
    }

    function range(r: NumRange) {
        return r && (r.max != null || r.min != null) ? `${r.min ?? ''}..${r.max ?? ''}` : null
    }
    function addMaybe(k: string, x: string | null) {
        if (x) {
            params.append(prefix + k, x)
        }
    }
    function add(k: string, x: string) {
        params.append(prefix + k, x)
    }

}

type UrlParseResult = {
    pairFamily: string | null
    pairType: string | null
    mode: 'basic' | 'ranges' | 'sql'
    filter: NucleotideFilterModel
    baselineFilter: NucleotideFilterModel | null | undefined
    stats: StatisticsSettingsModel | null
}

function parseRange(r: string | undefined | null): NumRange {
    if (!r) return {}
    const m = r.match(/^(-?\d+\.?\d*)?\.\.(-?\d+\.?\d*)?$/)
    if (!m) return {}
    const [_, min, max] = m
    return { min: min ? Number(min) : undefined, max: max ? Number(max) : undefined }
}

function parseRangeMaybe(r: string | undefined | null): NumRange | undefined {
    if (!r) return undefined
    const rr = parseRange(r)
    if (rr.min == null && rr.max == null)
        return undefined
    return rr
}

function trimArray(array: NumRange[]) {
    while (array.length && array.at(-1).max == null && array.at(-1).min == null) {
        array.pop()
    }
}

export function parseUrl(url: string): UrlParseResult {
    url = url.replace(/^[/]+/, '')
    const parts = url.split(/[/]+/)
    let pairType: string | null = null,
        pairFamily: string | null = null,
        m: RegExpMatchArray | null = null

    if (parts[0] && (m = parts[0].match(/^n?(?<fam>[ct][WHS][WHSB]a?)-(?<nt1>[ATUGC])-?(?<nt2>[ATGUC])$/i))) {
        pairFamily = m.groups?.fam
        pairType = `${m.groups?.fam}-${m.groups?.nt1}-${m.groups?.nt2}`
        parts.shift()
    }

    const f = new URLSearchParams(parts[0])
    const [filter, mode] = parseFilter(f)
    const [baselineFilter, _] = f.has('baseline_ds') || f.has('baseline_sql') ? parseFilter(f, undefined, 'baseline_') : [ undefined, "" ]

    const stats = parseStatsFromUrl(f)

    return { pairFamily, pairType, mode, filter, baselineFilter, stats }
}

function parseFilter(f: URLSearchParams, filter: NucleotideFilterModel | undefined = undefined, prefix=''): [NucleotideFilterModel, UrlParseResult['mode']] {
    filter ??= defaultFilter()
    filter.sql = f.get(`${prefix}sql`) ?? filter.sql
    if (f.has(`${prefix}ds`)) {
        let ds = f.get(`${prefix}ds`)
        if (/[DR]$/.test(ds)) {
            filter.dna = ds.endsWith('D')
            ds = ds.slice(0, -1)
        }
        filter.datasource = ds as NucleotideFilterModel['datasource']
    }
    else if (f.has(`${prefix}f`)) { // legacy links
        const fParam = f.get(`${prefix}f`)
        filter.filtered = fParam.includes('f')
        filter.includeNears = fParam.includes('n')
        filter.datasource = filter.filtered ? (filter.includeNears ? "fr3d-nf" : "fr3d-f") : (filter.includeNears ? "fr3d-n" : "fr3d")
        filter.dna = fParam.includes('D') ? true : fParam.includes('R') ? false : undefined
    }
    filter.bond_length = [{}, {}, {}, {}, {}, {}, {}, {}, {}, {}]
    filter.bond_donor_angle = [{}, {}, {}, {}, {}, {}, {}, {}, {}, {}]
    filter.bond_acceptor_angle = [{}, {}, {}, {}, {}, {}, {}, {}, {}, {}]
    filter.bond_plane_angle1 = [{}, {}, {}, {}, {}, {}, {}, {}, {}, {}]
    filter.bond_plane_angle2 = [{}, {}, {}, {}, {}, {}, {}, {}, {}, {}]
    for (let i = 0; i < 10; i++) {
        filter.bond_length[i] = parseRange(f.get(`${prefix}hb${i}_L`))
        filter.bond_donor_angle[i] = parseRange(f.get(`${prefix}hb${i}_DA`))
        filter.bond_acceptor_angle[i] = parseRange(f.get(`${prefix}hb${i}_AA`))
        filter.bond_plane_angle1[i] = parseRange(f.get(`${prefix}hb${i}_OOPA1`))
        filter.bond_plane_angle2[i] = parseRange(f.get(`${prefix}hb${i}_OOPA2`))
    }
    trimArray(filter.bond_length)
    trimArray(filter.bond_donor_angle)
    trimArray(filter.bond_acceptor_angle)
    trimArray(filter.bond_plane_angle1)
    trimArray(filter.bond_plane_angle2)

    for (const [k, v] of f.entries()) {
        if (k.startsWith(`${prefix}range_`)) {
            const column = k.slice(`${prefix}range_`.length)
            filter.other_column_range ??= {}
            filter.other_column_range[column] = parseRangeMaybe(v)
        }
    }

    if (f.has(`${prefix}condition`)) {
        filter.sql_conditions = f.getAll(`${prefix}condition`)
    }

    filter.coplanarity_angle = parseRangeMaybe(f.get(`${prefix}coplanar`) || f.get(`${prefix}coplanarity_a`))
    for (const prop of [ "coplanarity_edge_angle1", "coplanarity_edge_angle2", "coplanarity_shift1", "coplanarity_shift2", "yaw1", "yaw2", "pitch1", "pitch2", "roll1", "roll2" ]) {
        filter[prop] = parseRangeMaybe(f.get(`${prefix}${prop}`))
    }
    filter.resolution = parseRangeMaybe(f.get(`${prefix}resolution`) || f.get(`${prefix}resol`))
    filter.orderBy = f.get(`${prefix}order`) ?? filter.orderBy
    if (filter.orderBy) {
        // map expression back to orderBy ID
        filter.orderBy = orderByOptions.find(o => o.expr == filter.orderBy)?.id ?? filter.orderBy
    }
    const mode = filter.sql ? 'sql' : filter.bond_length.length || filter.bond_acceptor_angle.length || filter.bond_donor_angle.length || filter.coplanarity_angle || filter.bond_plane_angle1.length || filter.bond_plane_angle2.length || filter.coplanarity_edge_angle1 || filter.coplanarity_edge_angle2 || filter.coplanarity_shift1 || filter.coplanarity_shift2 || filter.pitch1 || filter.pitch2 || filter.yaw1 || filter.yaw2 || filter.roll1 || filter.roll2 || filter.resolution ? 'ranges' : 'basic'
    return [filter, mode]
}

function statPanelToStrings(params: URLSearchParams, ix: number, stat: StatPanelSettingsModel) {
    for (const [key, value] of Object.entries(statPresets)) {
        if (_.isEqual(value, stat)) {
            params.append(`st_${ix}`, "P" + key)
            return
        }
    }
    let prototypeVars = false
    for (const [key, value] of Object.entries(statPresets)) {
        if (_.isEqual(value.variables, stat.variables)) {
            params.append(`st_${ix}`, "P" + key)
            prototypeVars = true
        }
    }
    if (!prototypeVars) {
        params.append(`st_${ix}`, "T" + stat.type)
        if (stat.variables.some(v => v.filterSql)) {
            stat.variables.forEach(v => {
                params.append(`st_${ix}_v`, v.filterSql ? v.column + " WHERE " + v.filterSql : v.column)
            })
        } else {
            params.append(`st_${ix}_v`, stat.variables.map(v => v.column).join(','))
        }
        stat.variables.forEach((v, vix) => {
            if (v.label) {
                params.append(`st_${ix}_v${vix}_l`, v.label)
            }
        })
    }

    if (stat.title) {
        params.append(`st_${ix}_t`, stat.title)
    }
}

export function statsToUrl(params: URLSearchParams, stats: StatisticsSettingsModel) {
    if (!stats.enabled) {
        return
    }
    stats.panels.forEach((p, ix) => {
        statPanelToStrings(params, ix, p)
    })
}

export function parseStatsFromUrl(params: URLSearchParams): StatisticsSettingsModel | null {
    const keys = Array.from(params.keys()).filter(k => k.startsWith('st_'))
    if (keys.length == 0) {
        return null
    }
    const panelCount = Math.max(...keys.map(k => {
        const m = k.match(/^st_(\d+)/)
        return m ? Number(m[1])+1 : 0
    }))

    const result: StatisticsSettingsModel = { enabled: true, panels: [] }
    for (let i = 0; i < panelCount; i++) {
        const type = params.get(`st_${i}`)
        const panel = type.startsWith('P') ? _.cloneDeep(statPresets[type.slice(1)]) : { type: type.slice(1), variables: [], title:"" } as StatPanelSettingsModel
        const variables: VariableModel[] =
            params.getAll(`st_${i}_v`)
                .flatMap(v => v.includes(' WHERE ') ? [v] : v.split(','))
                .map(v => ({ label: "", column: v.split(' WHERE ')[0], filterSql: v.split(' WHERE ')[1] }))
        if (variables.length) {
            panel.variables = variables
        }
        panel.variables.forEach((v, vix) => {
            const label = params.get(`st_${i}_v${vix}_l`)
            if (label) {
                v.label = label
            }
        })
        const title = params.get(`st_${i}_t`)
        if (title) {
            panel.title = title
        }
        result.panels.push(panel)
    }
    return result
}

type UnwrapArray<T> = T extends Array<infer U> ? U : T

export const longBaseNames = { A: "adenine", T: "thymine", U: "uracil", G: "guanine", "C": "cytosine", DA: "adenine", DT: "thymine", DU: "uracil", DG: "guanine", DC: "cytosine" }
export const longNucleotideNames = {
    A: "adenosine", T: "thymidine", U: "uridine", G: "guanosine", "C": "cytidine",
    DA: "deoxy-adenosine", DT: "deoxy-thymidine", DU: "deoxy-uridine", DG: "deoxy-guanosine", DC: "deoxy-cytidine"
}

function capitalizeFirst(str: string | null | undefined) {
    return !str ? str : str.replace(/^./, x => x.toUpperCase())
}

function formatBaseNames(metadata: UnwrapArray<typeof metadataModule> | undefined, long: boolean, nucleotide: boolean = null) {
    if (!metadata) {
        return long ? [ 'base 1', 'base 2' ] : [ 'A', 'B' ]
    }
    nucleotide ??= long && metadata.atoms.some(a => a.some(a => a.endsWith("'") || a.startsWith("OP")))
    const [b1, b2] = metadata.pair_type[1].toUpperCase().split('-')
    if (long) {
        const dict = nucleotide ? longNucleotideNames : longBaseNames
        if (b1 == b2) {
            return [ `left ${dict[b1] ?? b1}`, `right ${dict[b2] ?? b2}` ]
        } else {
            return [ dict[b1] ?? b1, dict[b2] ?? b2 ]
        }
    } else {
        if (b1 == b2) {
            return [ `${b1}₁`, `${b2}₂` ]
        } else {
            return [ b1, b2 ]
        }
    }
}

function formatAtomNames(metadata: UnwrapArray<typeof metadataModule> | undefined, atoms: string[], long: boolean, comments?: (string | null)[]) {
    const [b1, b2] = formatBaseNames(metadata, long)
    let result = ''
    for (let i = 0; i < atoms.length; i++) {
        if (i > 0) {
            result += ', '
            if (long && i == atoms.length - 1)
                result += 'and '
        }
        if (/^[AB]/.test(atoms[i]) && atoms[i-1]?.[0] != atoms[i][0]) {
            result += (atoms[i][0] == 'A' ? b1 : b2) + ' '
        }

        result += atoms[i].replace(/^[AB]/, '')
        if (comments && comments[i])
            result += ' ' + comments[i]
    }
    return result
}

function baseNitrogen(base) {
    return base == 'G' || base == 'A' ? 'N9' : 'N1'
}

function formatEdge(e: string) {
    e = e.toUpperCase()
    return e == 'S' ? 'Sugar' :
           e == 'W' ? 'Watson-Crick' :
           e == 'H' ? 'Hoogsteen' :
           e == 'B' ? 'Bifurcated' :
           e
}

export function getColumnLabel(column: string, metadata: UnwrapArray<typeof metadataModule> | undefined, opt: { hideBondName?: boolean, hideParameterName?: boolean}={}): [string, string | null] {
    let m
    if ((m = /hb_(\d)_(\w+)/.exec(column))) {
        const [_, i, name] = m
        if (metadata?.labels[Number(i)]) {
            const hbondAtoms = metadata.atoms[i]
            const baseNames = formatBaseNames(metadata, true)
            const label = [
                opt.hideBondName ? null : metadata.labels[Number(i)],
                opt.hideParameterName ? null : ({
                    'length': "Length",
                    'donor_angle': "Donor-covalent ∡",
                    'acceptor_angle': "Acceptor-covalent ∡",
                    'donor_OOPA': 'Donor-plane ∡',
                    'acceptor_OOPA': 'Acceptor-plane ∡',
                    'OOPA1': `Left plane ∡`,
                    'OOPA2': `Right plane ∡`,
                }[name] ?? name)
            ].filter(x => x != null).join(' ')
            const tooltip =
                !hbondAtoms ? null :
                name == 'length' ? "Distance between the H-bond heavy atoms - " + formatAtomNames(metadata, hbondAtoms.slice(1, 3), true) :
                name == 'donor_angle' ? `Angle between heavy atoms ${formatAtomNames(metadata, hbondAtoms.slice(0, 3), true, [null, "(donor)", "(acceptor)"])}` :
                name == 'acceptor_angle' ? `Angle between heavy atoms ${formatAtomNames(metadata, hbondAtoms.slice(1, 4), true, ["(donor)", "(acceptor)", null])}` :
                // name == 'donor_OOPA' ? `Angle of ${metadata.labels[i]} and the ${baseNames[0]} plane` :
                // name == 'acceptor_OOPA' ? `Angle of ${metadata.labels[i]} and the ${baseNames[1]} (acceptor) plane` :
                name == 'OOPA1' ? `Angle of ${metadata.labels[i]} and the ${baseNames[0]} plane` :
                name == 'OOPA2' ? `Angle of ${metadata.labels[i]} and the ${baseNames[1]} plane` :
                null
            return [label, tooltip]
        }
    }
    if (column == "bogopropeller" || column == "coplanarity_angle") {
        return ["Coplanarity ∡", "Angle between planes of the pairing bases."]
    }
    if (column == 'coplanarity_shift1' || column == 'coplanarity_shift2') {
        let [b1, b2] = formatBaseNames(metadata, false, false)
        let [lb1, lb2] = formatBaseNames(metadata, true, false)
        if (column == 'coplanarity_shift1') {
            [b1, b2] = [b2, b1];
            [lb1, lb2] = [lb2, lb1]
        }
        return [`${b1} edge - ${b2} plane distance`, `Shortest distance between ${lb1} ${formatEdge(metadata.pair_type[0][column == 'coplanarity_shift2' ? 2 : 1])} edge and ${lb2} plane`]
    }
    if (column == 'coplanarity_edge_angle1' || column == 'coplanarity_edge_angle2') {
        const swap = column == 'coplanarity_edge_angle2'
        let [b1, b2] = formatBaseNames(metadata, false, false)
        let [lb1, lb2] = formatBaseNames(metadata, true, false)
        if (swap) {
            [b1, b2] = [b2, b1];
            [lb1, lb2] = [lb2, lb1]
        }
        return [`${b1} edge / ${b2} plane ∡`, `Angle between ${lb1} ${formatEdge(metadata.pair_type[0][swap ? 2 : 1])} edge and ${lb2} plane`]
    }
    if (column == "resolution") {
        return ["Structure resolution", null]
    }
    if (column == 'C1_C1_distance') {
        return ["C1'—C1' distance", null]
    }
    if (column == 'C1_C1_total_angle') {
        const [b1, b2] = formatBaseNames(metadata, true, true)
        const [n1, n2] = metadata.pair_type[1].toUpperCase().split('-').map(baseNitrogen)
        return [`C1'—${n1} / ${n2}—C1' ∡`, `Total angle between ${b1} C1'-${n1} bond and ${b2} ${n2}-C1' bond` ]
    }
    if ((m = /^C1_C1_(yaw|pitch|roll)(\d+)?$/.exec(column))) {
        const angleName = m[1]
        const [b1, b2] = formatBaseNames(metadata, true, false)
        const b = m[2] == '2' ? b2 : b1
        const [n1, n2] = metadata.pair_type[1].toUpperCase().split('-').map(baseNitrogen)
        return [`C1'—N ${capitalizeFirst(angleName)} ∡`, `${capitalizeFirst(angleName)} angle between C1'-${n1} bonds bond, relative to the ${b} nitrogen plane` ]
    }
    if (column == 'C1_C1_euler_phicospsi') {
        return [ `C1'—N Euler φ+cos(θ)ψ`, 'φ+cos(θ)ψ ≈ "total Z-axis rotation"' ]
    }
    if ((m = /^C1_C1_euler_(\w+)$/.exec(column))) {
        const angleName = m[1]
        const [b1, b2] = formatBaseNames(metadata, true, false)
        const [n1, n2] = metadata.pair_type[1].toUpperCase().split('-').map(baseNitrogen)
        const letter = { 'theta': 'θ', 'phi': 'φ', 'psi': 'ψ' }[angleName] ?? angleName
        return [`C1'—N Euler ${letter} ∡`, `${letter} (${angleName}) euler angle between C1'-${n1} bonds, relative to the ${b1} nitrogen plane` ]
    }
    if (column == 'rmsd_C1N_frames1' || column == 'rmsd_C1N_frames2') {
        let [b1, b2] = formatBaseNames(metadata, true, true)
        if (column == 'rmsd_C1N_frames2') {
            [b1, b2] = [b2, b1]
        }
        return [ `${column == 'rmsd_C1N_frames2' ? 'left' : 'right'} C1'—N RMSD`, `RMSD to exemplar basepair of ${b2} C1'—N reference frame atoms when aligned on ${b1} C1'—N` ]
    }
    if (column == 'rmsd_edge1' || column == 'rmsd_edge2') {
        let [b1, b2] = formatBaseNames(metadata, true, false)
        if (column == 'rmsd_edge2') {
            [b1, b2] = [b2, b1]
        }
        return [ `${column == 'rmsd_edge2' ? 'left' : 'right'} edge RMSD`, `RMSD to exemplar basepair of ${b2} ${formatEdge(metadata.pair_type[0][column == 'rmsd_edge2' ? 2 : 1])} edge, when aligned on ${b1}` ]
    }
    if (column == 'rmsd_all_base') {
        return [ 'All-atom RMSD', 'All-atom RMSD to exemplar basepair (only bases, not nucleotides)' ]
    }
    if (column == 'is_parallel') {
        return [ 'Strands are parallel', 'true if the basepair is between parallel strands, false if antiparallel, null if it cannot be determined' ]
    }
    if (column == 'is_dinucleotide') {
        return [ 'Is dinucleotide', 'true if the basepair are covalently linked consecutive nucleotides' ]
    }
    if (column == 'jirka_approves') {
        return [ "Passed quality filter", null ]
    }
    if (column == 'mode_deviations') {
        return [ "Sum of stddevs from KDE mode", "Sum of standard deviations from the peak of Kernel Density Estimate of selected parameters" ]
    }
    if (column == 'log_likelihood') {
        return [ "Log likelihood in KDE", "Sum of log likelihood in Kernel Density Estimate of selected parameters" ]
    }
    return null
}

export function hideColumn(column: string, metadata: UnwrapArray<typeof metadataModule> | undefined) {
    let m
    if (metadata && (m = /hb_(\d)_(\w+)/.exec(column))) {
        return Number(m[1]) >= metadata.atoms.length
    }
    return false
}

export function fillStatsLegends(stats: StatisticsSettingsModel, metadata: UnwrapArray<typeof metadataModule> | undefined): StatisticsSettingsModel {
    if (!stats.enabled || metadata == null) return stats

    const panels = stats.panels.map(p => fillStatLegends(p, metadata))
    return {...stats, panels }
}

export function fillStatLegends(p: HistogramSettingsModel, metadata: UnwrapArray<typeof metadataModule> | undefined): HistogramSettingsModel
export function fillStatLegends(p: KDE2DSettingsModel, metadata: UnwrapArray<typeof metadataModule> | undefined): KDE2DSettingsModel
export function fillStatLegends(p: StatPanelSettingsModel, metadata: UnwrapArray<typeof metadataModule> | undefined): StatPanelSettingsModel
export function fillStatLegends(p: StatPanelSettingsModel, metadata: UnwrapArray<typeof metadataModule> | undefined): StatPanelSettingsModel {
    const hideBondName = p.variables.length > 1 && p.variables.every(v => v.column.startsWith('hb_') && v.column.startsWith(p.variables[0].column.slice(0, 4)))
    const hideParameterName = p.variables.length > 1 && p.variables.every(v => v.column.startsWith('hb_') && v.column.endsWith(p.variables[0].column.slice(5)))
    const variables: VariableModel[] = p.variables.map(v => {
        if (!v.label) {
            const [label, tooltip] = getColumnLabel(v.column, metadata, { hideBondName, hideParameterName }) ?? [ v.label, v.tooltip ]
            v = { ...v, label, tooltip }
        }
        return v
    })
    if (p.type == 'kde2d') {
        return {...p, variables } as KDE2DSettingsModel
    } else if (p.type == 'histogram') {
        return {...p, variables } as HistogramSettingsModel
    }
    return {...(p as object), variables} as HistogramSettingsModel
}

export const statPresets: { [n: string]: StatPanelSettingsModel } = {
    "histL": { type: "histogram",
        title: "H-bond length (Å)",
        variables: [ { column: "hb_0_length", label: "" }, { column: "hb_1_length", label: "" }, { column: "hb_2_length", label: "" }, {column: "hb_3_length", label:""} ] },
    "histDA": { type: "histogram",
        title: "H-bond donor angle (°)",
        variables: [ { column: "hb_0_donor_angle", label: "" }, { column: "hb_1_donor_angle", label: "" }, { column: "hb_2_donor_angle", label: "" }, {column: "hb_3_donor_angle", label:""} ] },
    "histAA": { type: "histogram",
        title: "H-bond acceptor angle (°)",
        variables: [ { column: "hb_0_acceptor_angle", label: "" }, { column: "hb_1_acceptor_angle", label: "" }, { column: "hb_2_acceptor_angle", label: "" }, {column: "hb_3_acceptor_angle", label:""} ] }
}
