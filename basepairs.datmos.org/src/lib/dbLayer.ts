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
    dna?: true | false | undefined
    orderBy?: string
    filtered: boolean
    includeNears: boolean
}

export function defaultFilter(): NucleotideFilterModel {
    return { bond_acceptor_angle: [], bond_donor_angle: [], bond_length: [], filtered: true, includeNears: false }
}

function rangeToCondition(col: string, range: Range): string[] {
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
        if (filter.bond_length[i]) {
            conditions.push(...rangeToCondition(`hb_${i}_length`, filter.bond_length[i]))
        }
        if (filter.bond_donor_angle[i]) {
            conditions.push(...rangeToCondition(`hb_${i}_donor_angle`, filter.bond_donor_angle[i]))
        }
        if (filter.bond_acceptor_angle[i]) {
            conditions.push(...rangeToCondition(`hb_${i}_acceptor_angle`, filter.bond_acceptor_angle[i]))
        }
    }
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
}

function parseRange(r: string | undefined | null): Range {
    if (!r) return {}
    const m = r.match(/^(\d+)?\.\.(\d+)?$/)
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
    filter.orderBy = f.get(`order`)

    return { pairFamily, pairType, mode, filter }
}
