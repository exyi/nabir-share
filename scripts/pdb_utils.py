
import os, gzip, io, dataclasses
from typing import Any, Callable, Optional, TextIO, TypeVar, Union
import numpy as np


pdb_cache_dirs: list[str] = [ x for x in os.environ.get("PDB_CACHE_DIR", "").split(";") if x.strip() != "" ]
tmp_dir = "/tmp/pdb_files"

def _get_pdbid(file):
    pdbid = os.path.basename(file).split(".")[0]
    assert len(pdbid) == 4
    return pdbid

def _find_in_cache(pdb_id):
    for cache_dir in pdb_cache_dirs:
        if not os.path.isdir(cache_dir):
            continue
        all_files = os.listdir(cache_dir)
        all_files = { f.lower(): f for f in all_files }
        extensions = [ ".cif.gz", ".cif.zst", ".cif" ]

        # try direct file PDBID.cif
        for ext in extensions:
            if pdb_id.lower() + ext in all_files:
                return os.path.join(cache_dir, all_files[pdb_id.lower() + ext])

        # subdirectory PDBID[:2]/PDBID.cif
        subdir = all_files.get(pdb_id[:2].lower(), None)
        if subdir is not None:
            subdir = os.path.join(cache_dir, subdir)
            file = next((os.path.join(subdir, pdbfile) for pdbfile in os.listdir(subdir) if pdbfile.lower().startswith(pdb_id.lower() + ".")), None)
            if file is not None:
                return file
            
def _get_cache_write_path(pdb_id, compression = 'zst'):
    if len(pdb_cache_dirs) == 0:
        return None
    pdb_cache_dir = pdb_cache_dirs[-1]
    subdir = os.path.join(pdb_cache_dir, pdb_id[:2].lower())
    return os.path.join(subdir, pdb_id.lower() + ".cif." + compression)

def _find_or_add_cache(pdb_id):
    cache_file = _find_in_cache(pdb_id)
    if cache_file is not None:
        return cache_file
    else:
        write_path = _get_cache_write_path(pdb_id, compression='zst')
        if write_path is None:
            raise Exception("No cache directory")
        os.makedirs(os.path.dirname(write_path), exist_ok=True)
        with download_pdb(pdb_id) as f:
            write_path_ = write_path + ".part" + str(os.getpid())
            import zstandard
            with zstandard.open(write_path_, "wt", cctx=zstandard.ZstdCompressor(level=19)) as f2:
                f2.write(f.read())
            os.rename(write_path_, write_path)
        return write_path

def download_pdb(pdb_id: str) -> TextIO:
    import urllib.request
    # download from https://files.rcsb.org/download/XXXX.cif.gz
    url = f"https://files.rcsb.org/download/{pdb_id}.cif.gz"
    response = urllib.request.urlopen(url)
    return gzip.open(io.BytesIO(response.read()), "rt", encoding='utf-8-sig')

def get_pdb_file(file: Optional[str], pdb_id: Optional[str] = None) -> str:
    if file is not None:
        return file
    return _find_or_add_cache(pdb_id)
    
def open_pdb_file(file: Optional[str], pdb_id: Optional[str] = None) -> TextIO:
    if file is None and pdb_id is not None:
        if len(pdb_cache_dirs) == 0:
            return download_pdb(pdb_id)
        else:
            return open_pdb_file(_find_or_add_cache(pdb_id), pdb_id)

    elif isinstance(file, str):
        if file.endswith(".gz"):
            return gzip.open(file, "rt")
        elif file.endswith(".zst"):
            import zstandard
            return zstandard.open(file, "rt")
        elif file.endswith(".bz2"):
            import bz2
            return bz2.open(file, "rt")
        else:
            return open(file, "rt")
    elif file is None:
        raise Exception("No file specified")
    else:
        return file

def load_pdb_gemmi(file: Optional[str | TextIO], pdb_id: Optional[str] = None) -> 'Bio.PDB.Structure.Structure':
    import gemmi
    if isinstance(file, str) or file is None:
        with open_pdb_file(file, pdb_id) as f:
            return load_pdb_gemmi(f, pdb_id)
    else:
        raise NotImplementedError()
        # parser = Bio.PDB.MMCIFParser(QUIET=True)
        # structure = parser.get_structure(pdb_id, file)
        # return structure

def load_pdb(file: Optional[str | TextIO], pdb_id: Optional[str] = None):
    import Bio.PDB.Structure
    if isinstance(file, str) or file is None:
        with open_pdb_file(file, pdb_id) as f:
            return load_pdb(f, pdb_id)
    else:
        parser = Bio.PDB.MMCIFParser(QUIET=True)
        structure = parser.get_structure(pdb_id, file)
        return structure

@dataclasses.dataclass
class SymmetryOperation:
    id: str
    pdbname: str
    triplet: str
    rotation: np.ndarray
    translation: np.ndarray


@dataclasses.dataclass
class StructureData:
    assembly: list[SymmetryOperation]
    organism: Optional[str] = None
    organism_id: Optional[int] = None

def load_sym_data(file: Optional[str | TextIO], pdb_id: Optional[str] = None) -> StructureData:
    from mmcif.io.PdbxReader import PdbxReader
    if isinstance(file, str) or file is None:
        with open_pdb_file(file, pdb_id) as f:
            return load_sym_data(f, pdb_id)
    else:
        datablocks = []
        PdbxReader(file).read(datablocks)
        result = StructureData([])

        assembly_obj = datablocks[0].getObj("pdbx_struct_oper_list")
        if assembly_obj:
            idx = assembly_obj.getAttributeIndexDict()
            for row in assembly_obj.data:
                result.assembly.append(SymmetryOperation(
                    id=str(row[idx["id"]]),
                    pdbname=row[idx["name"]],
                    triplet=row[idx["symmetry_operation"]],
                    rotation=np.array([
                        [ float(row[idx[f"matrix[1][1]"]]), float(row[idx[f"matrix[1][2]"]]), float(row[idx[f"matrix[1][3]"]]) ],
                        [ float(row[idx[f"matrix[2][1]"]]), float(row[idx[f"matrix[2][2]"]]), float(row[idx[f"matrix[2][3]"]]) ],
                        [ float(row[idx[f"matrix[3][1]"]]), float(row[idx[f"matrix[3][2]"]]), float(row[idx[f"matrix[3][3]"]]) ],
                    ]),
                    translation=np.array([ float(row[idx[f"vector[1]"]]), float(row[idx[f"vector[2]"]]), float(row[idx[f"vector[3]"]]) ])
                ))

        if obj := datablocks[0].getObj("pdbx_entity_src_syn"):
            idx = obj.getAttributeIndexDict()
            result.organism = obj.data[0][idx["organism_scientific"]]
            result.organism_id = int(obj.data[0][idx["pdbx_src_id"]]) if obj.data[0][idx["pdbx_src_id"]] else None
        return result

