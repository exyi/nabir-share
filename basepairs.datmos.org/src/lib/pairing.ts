export type NucleotideId = {
    pdbid: string;
    model: number;
    chain: string;
    resname?: string;
    resnum: number;
    altloc: string;
    inscode: string;
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
}

export type PairingInfo = {
    id: PairId
    hbonds?: HydrogenBondInfo[]
    coplanarity?: number
    originalRow?: any
}

export function normalizePairFamily(f: string) {
    return f.replace('s', 'S').replace('w', 'W').replace('h', 'H') as PairingFamily
}

export function tryParsePairingType(s: string | [ PairingFamily, string ]): [ PairingFamily, string ] | undefined {
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
