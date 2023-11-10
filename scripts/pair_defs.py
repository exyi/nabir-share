import itertools
import os, math, re, csv, dataclasses
from typing import Union

@dataclasses.dataclass(frozen=True)
class PairType:
    type: str
    bases: tuple[str, str]

    def __post_init__(self):
        # assert len(self.type) == 3
        # assert len(self.bases) == 2
        # assert self.type[0] in "ct"
        # assert self.type[1] in "WSH"
        # assert self.type[2] in "WSH"
        # assert self.bases[0] in "ACGU"
        # assert self.bases[1] in "ACGU"
        if len(self.type) != 3 or len(self.bases) != 2 or self.type[0] not in "ct" or self.type[1] not in "WSH" or self.type[2] not in "WSH" or (len(self.bases[0]) == 1 and self.bases[0] not in "ACGU") or (len(self.bases[1]) == 1 and self.bases[1] not in "ACGU"):
            raise ValueError(f"Invalid pair type: PairType({repr(self.type)}, {self.bases})")
    def __str__(self) -> str:
        return f"{self.type}-{'-'.join(self.bases)}"
    def __repr__(self) -> str:
        return f"PairType({repr(self.type)}, {self.bases})"
    def to_tuple(self) -> tuple[str, str]:
        return (self.type, self.bases[0] + "-" + self.bases[1])
    def swap(self) -> 'PairType':
        return PairType(self.type[0] + self.type[2] + self.type[1], (self.bases[1], self.bases[0]))
    def is_preferred_orientation(self) -> bool:
        return is_preferred_pair_type_orientation(self.to_tuple())
    def is_swappable(self):
        if self.to_tuple() in hbonding_atoms and self.swap().to_tuple() in hbonding_atoms:
            # both definitions exist,
            return False
        return True
    @staticmethod
    def from_tuple(t: tuple[str, str]) -> "PairType":
        return PairType(t[0], tuple(t[1].split("-"))) # type: ignore
    
    @staticmethod
    def parse(s: str) -> "PairType":
        raise NotImplementedError("PairType.parse is not implemented")

def read_pair_definitions(file = os.path.join(os.path.dirname(__file__), "H_bonding_Atoms_from_Isostericity_Table.csv")) -> dict[tuple[str, str], list[tuple[str, str, str, str]]]:
    with open(file, "r") as f:
        reader = csv.reader(f)
        lines = list(reader)
    header_i = next(i for i, line in enumerate(lines) if len(line) > 3 and line[0] == 'TYPE')
    header = lines[header_i]
    data = lines[header_i + 1:]
    def translate_pair_type(line):
        m = re.match(r"(cis |trans )([WSH])/([WSH])", line[0].strip(), re.IGNORECASE)
        if not m:
            # print("WARNING: Invalid pair type: " + line[0])
            return None
        ct = { "c": "c", "t": "t", "cis":"c", "trans": "t" }[m.group(1)[0].strip().lower()]
        return f"{ct}{m.group(2).upper()}{m.group(3).upper()}"
    
    result_mapping = {}
    for line in data:
        res1, res2 = line[1].strip(), line[2].strip()
        pt = (translate_pair_type(line), res1 +"-"+ res2)
        if pt[0] is None:
            continue
        if pt not in result_mapping:
            result_mapping[pt] = []
        assert line[5] == '---', line
        atom1, atom2, atom3, atom4 = line[3].strip(), line[4].strip(), line[6].strip(), line[7].strip()
        if atom1 == '-':
            # atom2 is the acceptor, atom3 is hydrogen and atom4 is the donor
            b = ("B"+get_angle_ref_atom(res2, atom4), "B"+atom4, "A"+atom2, "A"+get_angle_ref_atom(res1, atom2))
        else:
            assert atom4 == '-'
            b = ("A"+get_angle_ref_atom(res1, atom1), "A"+atom1, "B"+atom3, "B"+get_angle_ref_atom(res2, atom3))
        result_mapping[pt].append(b)
    return result_mapping

def get_angle_ref_atom(res: str, atom: str)-> str:
    res = {"U":"T"}.get(res, res)
    neighbors = set(itertools.chain(*[ edge for edge in atom_connectivity[res] if atom in edge ]))
    neighbors.remove(atom)
    assert len(neighbors) >= 1, f"No neighbors found in {res} for atom {atom}"
    return max(neighbors)

def hbond_swap_nucleotides(hbond: tuple[str, str, str, str]) -> tuple[str, str, str, str]:
    return tuple(
        atom.replace("A", "ⒷⒷⒷⒷⒷⒷ").replace("B", "A").replace("ⒷⒷⒷⒷⒷⒷ", "B")
        for atom in hbond
    ) # type: ignore

def swap_pair_type(pair_type: tuple[str, str]) -> tuple[str, str]:
    t, bases = pair_type
    bases = bases.split("-")
    assert len(pair_type) == 2 and len(t) == 3 and len(bases) == 2
    assert t[0] in "ct"
    return (t[0] + t[2] + t[1], bases[1] + "-" + bases[0])

def is_preferred_pair_type_orientation(pair_type: tuple[str, str]) -> bool:
    other = swap_pair_type(pair_type)
    if other == pair_type:
        return True
    if other[0] == pair_type[0]:
        # symetric pair type like cWW, cHH, ...
        # we prefer A > G > C > U (who knows why)
        preference = [ "A", "G", "C", "U" ]
        bases = pair_type[1].split("-")
        if bases[0] not in preference and bases[1] not in preference:
            # weird bases, alphabetical order
            return bases[0] <= bases[1]
        elif bases[0] not in preference:
            return False
        elif bases[1] not in preference:
            return True

        return preference.index(pair_type[1][0]) <= preference.index(other[1][0])
    else:
        # non-symetric pair type like cWH, cWS, ...
        return pair_type[0] in pair_types

# all pair type ordered according to The Paper
pair_types = [
    "cWW", "tWW",
    "cWH", "tWH",
    "cWS", "tWS",
    "cHH", "tHH",
    "cHS", "tHS",
    "cSS", "tSS"
]

sugar_atom_connectivity = [
    ("P", "OP2"), ("P", "OP1"), ("P", "O5'"),
    ("O5'", "C5'"), ("C5'", "C4'"),
    ("C4'", "O4'"), ("C4'", "C3'"),
    ("C3'", "O3'"), ("C3'", "C2'"),
    ("C2'", "O2'"), ("C2'", "C1'"),
    ("C1'", "O4'")
]
purine_atom_connectivity = [
    *sugar_atom_connectivity,
    ("C1'", "N9"), ("N9", "C8"), ("C8", "N7"),
    ("N9", "C4"), ("N7", "C5"),
    ("C4", "C5"), ("C5", "C6"), ("C6", "N1"), ("N1", "C2"), ("C2", "N3"), ("N3", "C4")
]
pyrimidine_atom_connectivity = [
    *sugar_atom_connectivity,
    ("C1'", "N1"), ("N1", "C2"), ("C2", "N3"),
    ("N3", "C4"), ("C4", "C5"), ("C5", "C6"),
    ("N1", "C6"),
]
atom_connectivity = {
    "A": [
        *purine_atom_connectivity,
        ("C6", "N6")
    ],
    "G": [
        *purine_atom_connectivity,
        ("C6", "O6"), ("C2", "N2")
    ],
    "C": [
        *pyrimidine_atom_connectivity,
        ("C2", "O2"),
        ("C4", "N4"),
        ("C5", "C7") # 5-methyl cytosine
    ],
    "T": [
        *pyrimidine_atom_connectivity,
        ("C2", "O2"),
        ("C4", "O4"),
        ("C5", "C7")
    ]
}

my_hbonding_atoms: dict[tuple[str, str], list[tuple[str, str, str, str]]] = {
    ('cWW', 'C-G'): [
        ('AC4', 'AN4', 'BO6', 'BC6'),
        ('BC6', 'BN1', 'AN3', 'AC2'),
        ('BC2', 'BN2', 'AO2', 'AC2'),
    ],
    ('cWW', 'A-U'): [
        ('AC6', 'AN6', 'BO4', 'BC4'),
        ('BC2', 'BN3', 'AN1', 'AC2') # TODO: change BC2 to BC4 and AC2 to AC6
    ],
    ('cWW', 'G-U'): [
        ('BC4', 'BN3', 'AO6', 'AC6'),
        ('AC2', 'AN1', 'BO2', "BC2"),
        ('AC2', 'AN2', "BO2'", "BC2'"), # ??? pres vodu
    ],
    ('tHS', 'A-G'): [
        ('BC2', 'BN2', 'AN7', 'AC8'),
        ('AC6', 'AN6', 'BN3', 'BC4'),
        ('AC6', 'AN6', "BO2'", "BC2'"),
    ],
    # ('cSS', 'A-C'): [ ], # kokotina
    ('tSS', 'A-G'): [
        ("BC2", "BN2", "AN3", "AC2"),
        ("AN3", "AC2", "BN3", "BC4"),
        ("BC2'", "BO2'", "AN1", "AC6")
    ],
    ('tHW', 'A-U'): [
        ("BC4", "BN3", "AN7", "AC8"),
        ("AC6", "AN6", "BO2", "BC2"),
    ],
    ('cWW', 'U-U'): [
        ("BC4", "BN3", "AO4", "AC4"),
        ("AC4", "AN3", "BO2", "BC2"),
        ("BC2'", "BO2'", "AO2", "AC2"), # pres vodu
    ],
    ('tHH', 'A-A'): [
        ("AC6", "AN6", "BN7", "BC8"),
        ("BC6", "BN6", "AN7", "AC8"),
        ("BC6", "BN6", "AOP2", "AP"),
    ],
    ('cSS', 'A-G'): [
        ("BC2", "BN2", "AN3", "AC4"),
        ("AC2'", "AO2'", "BN3", "BC4"),
        ("BC2'", "BO2'", "AO2'", "AC2'"),
    ],
    ('cWW', 'A-G'): [
        ('AC6', 'AN6', 'BO6', 'BC6'),
        ('BC6', 'BN1', 'AN1', 'AC6'),
        # ('BC2', 'BN2', 'AC2', 'AN3'),
    ],
    ('tWW', 'A-A'): [
        ("AC6", "AN6", "BN1", "BC6"),
        ("BC6", "BN6", "AN1", "AC6")
    ],
    # ('tSS', 'A-C'): [ ], # kokotina
    ('tWS', 'A-G'): [
        ("BC2", "BN2", "AN1", "AC6"),
        ("AC6", "AN6", "BN3", "BC4"),
        ("AC6", "AN6", "BO2'", "BC2'"),
    ],
    # ('cSS', 'A-U'): [ ], # kokotina
    # ('cSH', 'G-U'): [ ], # kokotina
    ('cSS', 'A-A'): [
        ("BN3", "BC2", "AN3", "AC4"),
        ("AC2'", "AO2'", "BN3", "BC4"),
        ("BC2'", "BO2'", "AO2'", "AC2'"),
    ],
    ('tSS', 'G-G'): [
        ("BC2", "BN2", "AN3", "AC4"),
        ("AC2", "AN2", "BN3", "BC4"),
        ("AC2", "AN2", "BO2'", "BC2'"),
        # ("BC2", "BN2", "AO2'", "AC2'"),
    ],
    ('cWH', 'G-G'): [
        ("AC2", "AN2", "BN7", "BC8"),
        ("AC6", "AN1", "BO6", "BC6")
    ],
    ('cWW', 'C-U'): [
        ("AC4", "AN4", "BO4", "BC4"),
        ("BC4", "BN3", "AN3", "AC4"),
    ],
    # ('cWH', 'C-U'): [
    #     ("AC4", "AN4", "BO4", "BC4"),
    #     ("BC2", "BN3", "AN3", "AC4"),
    # ],
}

# hbonding_atoms = my_hbonding_atoms
hbonding_atoms = read_pair_definitions()


def is_bond_hidden(pair_type, b) -> bool:
    return b[1][1] == 'C' or b[2][1] == 'C' or \
        b[1].endswith("'") or b[2].endswith("'")

def get_hbonds(pair_type: Union[tuple[str, str], PairType], throw=True) -> list[tuple[str, str, str, str]]:
    if isinstance(pair_type, PairType):
        pair_type = pair_type.to_tuple()
    if pair_type in hbonding_atoms:
        return hbonding_atoms[pair_type]
    
    pair_type_s = swap_pair_type(pair_type)
    if pair_type_s in hbonding_atoms:
        return [ hbond_swap_nucleotides(hbond) for hbond in hbonding_atoms[pair_type_s] ]
    
    if throw:
        raise KeyError(f"Pair type {pair_type} not found in hbonding_atoms")
    else:
        return []

_resname_map = {
    'DT': 'U',
    'DC': 'C',
    'DA': 'A',
    'DG': 'G',
    'DU': 'U',
    'T': 'U',
}
def map_resname(resname: str) -> str:
    return _resname_map.get(resname.upper(), resname)

