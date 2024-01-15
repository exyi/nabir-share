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
