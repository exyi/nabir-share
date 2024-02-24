import type metadataModule from './metadata'
import _ from 'lodash'


export type Range = {
    min?: number
    max?: number
}
export type NucleotideFilterModel = {
    sql?: string
    bond_length: (Range)[]
    bond_donor_angle: (Range)[]
    bond_acceptor_angle: (Range)[]
    coplanarity?: Range
    resolution?: Range
    dna?: true | false | undefined
    orderBy?: string
    filtered: boolean
    includeNears: boolean
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

export function defaultFilter(): NucleotideFilterModel {
    return { bond_acceptor_angle: [], bond_donor_angle: [], bond_length: [], filtered: true, includeNears: false }
}

function rangeToCondition(col: string, range: Range | undefined | null): string[] {
    if (range == null) return []
    const r = []
    if (range.min != null) {
        r.push(`${col} >= ${range.min}`)
    }
    if (range.max != null) {
        r.push(`${col} <= ${range.max}`)
    }
    return r
}

export function filterToSqlCondition(filter: NucleotideFilterModel) {
    const conditions = []
    if (filter.filtered) {
        conditions.push(`jirka_approves`)
    }
    for (let i = 0; i < 3; i++) {
        conditions.push(...rangeToCondition(`hb_${i}_length`, filter.bond_length[i]))
        conditions.push(...rangeToCondition(`hb_${i}_donor_angle`, filter.bond_donor_angle[i]))
        conditions.push(...rangeToCondition(`hb_${i}_acceptor_angle`, filter.bond_acceptor_angle[i]))
    }
    conditions.push(...rangeToCondition(`resolution`, filter.resolution))
    if (filter.coplanarity) {
        conditions.push(...rangeToCondition(`bogopropeller`, filter.coplanarity))
    }
    if (filter.dna != null) {
        if (filter.dna)
            conditions.push(`(res1 LIKE 'D%' OR res2 LIKE 'D%')`)
        else
            conditions.push(`(res1 NOT LIKE 'D%' OR res2 NOT LIKE 'D%')`)
    }
    return conditions
}

export function makeSqlQuery(filter: NucleotideFilterModel, from: string, limit?: number) {
    // if (filter.sql) {
    //     return filter.sql
    // }
    const conditions = filterToSqlCondition(filter)
    const where = conditions.map(c => /\b(select|or)\b/.test(c) ? `(${c})` : c).join(' AND ')

    let query = `SELECT * FROM ${from}`
    if (where) {
        query += ` WHERE ${where}`
    }
    if (filter.orderBy) {
        query += ` ORDER BY ${filter.orderBy}`
    }
    if (limit) {
        query += ` LIMIT ${limit}`
    }
    return query
}

export function aggregateTypesQuery(query: string, type = "type", res1 = "res1", res2 = "res2") {
    return `
        SELECT concat(${type}, '-', ltrim(${res1}, 'D'), '-', ltrim(${res2}, 'D')) as type,
               COUNT(*) AS count
        FROM (${query})
        GROUP BY ${type}, ltrim(${res1}, 'D'), ltrim(${res2}, 'D')
        ORDER BY COUNT(*) DESC`
}

export function aggregatePdbCountQuery(query: string, pdbid = "pdbid") {
    return `
        SELECT ${pdbid} as pdbid, COUNT(*) AS count
        FROM (${query})
        GROUP BY ${pdbid}
        ORDER BY COUNT(*) DESC`
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


export function filterToUrl(filter: NucleotideFilterModel, mode = 'ranges') {
    if (mode == 'sql') {
        return new URLSearchParams({ sql: filter.sql })
    }

    function range(r: Range) {
        return r && (r.max || r.min) ? `${r.min ?? ''}..${r.max ?? ''}` : null
    }
    function addMaybe(k, x: string | null) {
        if (x) {
            params.append(k, x)
        }
    }
    const params = new URLSearchParams()
    if (!filter.filtered || filter.includeNears || filter.dna != null) {
        params.set('f', (filter.filtered ? 'f' : '') + (filter.includeNears ? 'n' : '') + (filter.dna == true ? 'D' : filter.dna == false ? 'R' : ''))
    }
    for (let i = 0; i < 3; i++) {
        addMaybe(`hb${i}_L`, range(filter.bond_length[i]))
        addMaybe(`hb${i}_DA`, range(filter.bond_donor_angle[i]))
        addMaybe(`hb${i}_AA`, range(filter.bond_acceptor_angle[i]))
    }
    addMaybe(`coplanar`, range(filter.coplanarity))
    addMaybe(`resol`, range(filter.resolution))
    if (filter.orderBy) {
        addMaybe(`order`, filter.orderBy)
    }
    return params
}

type UrlParseResult = {
    pairFamily: string | null
    pairType: string | null
    mode: 'ranges' | 'sql'
    filter: NucleotideFilterModel
    stats: StatisticsSettingsModel | null
}

function parseRange(r: string | undefined | null): Range {
    if (!r) return {}
    const m = r.match(/^(-?\d+\.?\d*)?\.\.(-?\d+\.?\d*)?$/)
    if (!m) return {}
    const [_, min, max] = m
    return { min: min ? Number(min) : undefined, max: max ? Number(max) : undefined }
}

function trimArray(array: Range[]) {
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

    const filter = defaultFilter()
    const f = new URLSearchParams(parts[0])
    filter.sql = f.get('sql')
    const mode = filter.sql ? 'sql' : 'ranges'
    if (f.has('f')) {
        filter.filtered = f.get('f').includes('f')
        filter.includeNears = f.get('f').includes('n')
        filter.dna = f.get('f').includes('D') ? true : f.get('f').includes('R') ? false : undefined
    }
    filter.bond_length = [{}, {}, {}, {}, {}, {}, {}, {}, {}, {}]
    filter.bond_donor_angle = [{}, {}, {}, {}, {}, {}, {}, {}, {}, {}]
    filter.bond_acceptor_angle = [{}, {}, {}, {}, {}, {}, {}, {}, {}, {}]
    for (let i = 0; i < 10; i++) {
        filter.bond_length[i] = parseRange(f.get(`hb${i}_L`))
        filter.bond_donor_angle[i] = parseRange(f.get(`hb${i}_DA`))
        filter.bond_acceptor_angle[i] = parseRange(f.get(`hb${i}_AA`))
    }
    trimArray(filter.bond_length)
    trimArray(filter.bond_donor_angle)
    trimArray(filter.bond_acceptor_angle)

    filter.coplanarity = parseRange(f.get(`coplanar`))
    filter.resolution = parseRange(f.get(`resolution`) || f.get(`resol`))
    filter.orderBy = f.get(`order`)

    const stats = parseStatsFromUrl(f)

    return { pairFamily, pairType, mode, filter, stats }
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
            stat.variables.forEach((v, vix) => {
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
                }[name] ?? name)
            ].filter(x => x != null).join(' ')
            const tooltip =
                !hbondAtoms ? null :
                name == 'length' ? "Distance between the H-bond heavy atoms - " + formatAtomNames(metadata, hbondAtoms.slice(1, 3), true) :
                name == 'donor_angle' ? `Angle between heavy atoms ${formatAtomNames(metadata, hbondAtoms.slice(0, 3), true, [null, "(donor)", "(acceptor)"])}` :
                name == 'acceptor_angle' ? `Angle between heavy atoms ${formatAtomNames(metadata, hbondAtoms.slice(1, 4), true, ["(donor)", "(acceptor)", null])}` :
                name == 'donor_OOPA' ? `Angle of ${metadata.labels[i]} and the ${baseNames[0]} (donor) plane` :
                name == 'acceptor_OOPA' ? `Angle of ${metadata.labels[i]} and the ${baseNames[1]} (acceptor) plane` :
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
        return ["C1'-C1' distance", null]
    }
    if (column == 'C1_C1_total_angle') {
        const [b1, b2] = formatBaseNames(metadata, true, true)
        const [n1, n2] = metadata.pair_type[1].toUpperCase().split('-').map(baseNitrogen)
        return [`C1'-${n1} / ${n2}-C1' ∡`, `Total angle between ${b1} C1'-${n1} bond and ${b2} ${n2}-C1' bond` ]
    }
    if ((m = /^C1_C1_(yaw|pitch|roll)(\d+)?$/.exec(column))) {
        const angleName = m[1]
        const [b1, b2] = formatBaseNames(metadata, true, false)
        const b = m[2] == '2' ? b2 : b1
        const [n1, n2] = metadata.pair_type[1].toUpperCase().split('-').map(baseNitrogen)
        return [`C1'-${n1} / ${n2}-C1' ${capitalizeFirst(angleName)} ∡`, `${capitalizeFirst(angleName)} angle between C1'-${n1} bonds bond, relative to the ${b} nitrogen plane` ]
    }
    if ((m = /^C1_C1_euler_(\w+)$/.exec(column))) {
        const angleName = m[1]
        const [b1, b2] = formatBaseNames(metadata, true, false)
        const [n1, n2] = metadata.pair_type[1].toUpperCase().split('-').map(baseNitrogen)
        const letter = { 'theta': 'θ', 'phi': 'φ', 'psi': 'ψ' }[angleName] ?? angleName
        return [`C1'-${n1} / ${n2}-C1' Euler ${letter} ∡`, `${letter} (${angleName}) euler angle between C1'-${n1} bonds, relative to the ${b1} nitrogen plane` ]
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
