export type NucleotideId = {
    pdbid: string;
    model: number;
    chain: string;
    resname?: string;
    resnum: number;
    altloc: string;
    inscode: string;
}
export type PairingType = `${'t' | 'c'}${'W' | 'H' | 'S'}${'W' | 'H' | 'S'}`
export type PairId = {
    nt1: NucleotideId;
    nt2: NucleotideId;
    pairingType?: [PairingType, string]
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
}


export function parsePairingType(s): [ PairingType, string ] {
    if (Array.isArray(s))
        return s as [ PairingType, string ]

    const m = s.match(/^(?<bases>[ACTG]-[ACTG])-(?<pt>[tc][WHS][WHS])$/)
    if (m) {
        return [ m.groups.pt as PairingType, m.groups.bases ]
    }
    throw new Error(`Invalid pairing type: ${s}`)
}
