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
}

function rangeToCondition(col, range: Range): string[] {
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

export function makeSqlQuery(filter: NucleotideFilterModel, from: string) {
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
    return query
}
