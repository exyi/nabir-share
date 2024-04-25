import metadata from "./metadata"
import type * as arrow from 'apache-arrow'

export type NucleotideId = {
    pdbid: string
    model: number
    chain: string
    resname?: string
    resnum: number
    altloc: string
    inscode: string
    symop?: string
}
export type PairingFamily = `${'n' | ''}${'t' | 'c'}${'W' | 'H' | 'S'}${'W' | 'H' | 'S'}${'a' | ''}`
export type PairId = {
    nt1: NucleotideId;
    nt2: NucleotideId;
    pairingType?: [PairingFamily, string]
}

export type HydrogenBondInfo = {
    label?: string
    length: number
    acceptorAngle: number
    donorAngle: number
    OOPA1: number
    OOPA2: number
}

export type PairingInfo = {
    id: PairId
    hbonds?: HydrogenBondInfo[]
    coplanarity?: number
    comparison?: boolean
    originalRow?: any
}

export function normalizePairType(t: string) {
    const [family, bases] = parsePairingType(t)
    return `${normalizePairFamily(family)}-${bases.toUpperCase()}`
}

export function normalizePairFamily(f: string) {
    return f.replace('s', 'S').replace('w', 'W').replace('h', 'H') as PairingFamily
}

export function tryParsePairingType(s: string | [ PairingFamily, string ]): [ PairingFamily, string ] | undefined {
    if (!s)
        return undefined
    if (Array.isArray(s))
        return s as [ PairingFamily, string ]

    const m =
        s.match(/^(?<bases>[ACTUG]-[ACTUG])-(?<n>n?)(?<pt>[tc][WHS][WHSwhs]a?)$/) ||
        s.match(/^(?<n>n?)(?<pt>[tc][WHS][WHSwhs]a?)-(?<bases>[ACTUG]-[ACTUG])$/)
    if (m) {
        if (m.groups == null) {
            return null
        }
        return [ normalizePairFamily(m.groups.pt), m.groups.bases ]
    }
    return null
}

export function parsePairingType(s: string | [ PairingFamily, string ]): [ PairingFamily, string ] {
    const result = tryParsePairingType(s)
    if (result == null) {
        throw new Error(`Invalid pairing type: ${s}`)
    }
    return result
}



export function compareFamilies(a: string, b: string) {
    function score(edge: string) {
        switch (edge.toUpperCase()) {
            case 'W': return 0;
            case 'H': return 1;
            case 'S': return 2;
            default: return 3;
        }
    }
    function edgeCompare(a: string, b: string) {
        return score(a) - score(b) || a.toUpperCase().localeCompare(b.toUpperCase())
    }
    function trimN(x: string) {
        if (x[0] == 'n') {
            return x.substr(1)
        }
        return x
    }

    const an = trimN(a)
    const bn = trimN(b)

    return edgeCompare(an[1], bn[1]) || edgeCompare(an[2], bn[2]) || an[0].localeCompare(bn[0]) || Number(a[0] != 'n') - Number(b[0] != 'n')
}
export function getMetadata(pairType) {
    if (pairType == null) return undefined
    if (Array.isArray(pairType)) {
        pairType = pairType.join('-')
    }
    pairType = pairType.toLowerCase()
    return metadata.find(m => m.pair_type.join('-').toLowerCase() == pairType)
}


export function* convertQueryResults(rs: arrow.Table, pairingType, limit=undefined): Generator<PairingInfo> {
    function convName(name) {
      if (name == null) return undefined
      name = String(name).trim()
      if (name == '' || name == '?' || name == '\0') return undefined
      return name
    }

    const meta = getMetadata(pairingType)

    let count = 0
    for (const r of rs) {
      if (limit != null && count >= limit)
        break
      const pdbid = r.pdbid, model = Number(r.model ?? 1)
      const nt1: NucleotideId = { pdbid, model, chain: convName(r.chain1), resnum: Number(r.nr1), resname: convName(r.res1), altloc: convName(r.alt1), inscode: convName(r.ins1), symop: r.symmetry_operation1 }
      const nt2: NucleotideId = { pdbid, model, chain: convName(r.chain2), resnum: Number(r.nr2), resname: convName(r.res2), altloc: convName(r.alt2), inscode: convName(r.ins2), symop: r.symmetry_operation2 }

      const id: PairId = { nt1, nt2, pairingType }
      let hbonds: HydrogenBondInfo[] | undefined
      if ("hb_0_length" in r) {
        hbonds = []
        for (let i = 0; i <= 4; i++) {
          const length = r[`hb_${i}_length`]
          if (length == null) break
          // eslint-disable-next-line no-inner-declarations
          function numberMaybe(x) { return x ? Number(x) : null }
          hbonds.push({
            length: numberMaybe(length),
            donorAngle: numberMaybe(r[`hb_${i}_donor_angle`]),
            acceptorAngle: numberMaybe(r[`hb_${i}_acceptor_angle`]),
            OOPA1: numberMaybe(r[`hb_${i}_OOPA1`]),
            OOPA2: numberMaybe(r[`hb_${i}_OOPA2`]),
            label: meta?.labels[i],
          })
        }
      }
      const coplanarity = r.bogopropeller ?? r.coplanarity_angle
      const { comparison_in_baseline, comparison_in_current } = r
      const comparison = (comparison_in_baseline || false) != (comparison_in_current || false) ? Boolean(comparison_in_current || false) : undefined
      yield { id, hbonds, coplanarity, comparison, originalRow: r }
      count++
    }
  }
