
import os, gzip, io, dataclasses
from typing import Any, Callable, Optional, TextIO, TypeVar, Union
import numpy as np


pdb_cache_dirs: list[str] = [ x for x in os.environ.get("PDB_CACHE_DIR", "").split(";") if x.strip() != "" ]
tmp_dir = "/tmp/pdb_files"

def _get_pdbid(file):
    """Extract PDBID from file name."""
    pdbid = os.path.basename(file).split(".")[0]
    assert len(pdbid) == 4
    return pdbid

def _find_in_cache(pdbid):
    """Try to find cif file in the `pdb_cache_dirs`, returns the file path"""
    for cache_dir in pdb_cache_dirs:
        if not os.path.isdir(cache_dir):
            continue
        all_files = os.listdir(cache_dir)
        all_files = { f.lower(): f for f in all_files }
        extensions = [ ".cif.gz", ".cif.zst", ".cif" ]

        # try direct file PDBID.cif
        for ext in extensions:
            if pdbid.lower() + ext in all_files:
                return os.path.join(cache_dir, all_files[pdbid.lower() + ext])

        # subdirectory PDBID[:2]/PDBID.cif
        subdir = all_files.get(pdbid[:2].lower(), None)
        if subdir is not None:
            subdir = os.path.join(cache_dir, subdir)
            file = next((os.path.join(subdir, pdbfile) for pdbfile in os.listdir(subdir) if pdbfile.lower().startswith(pdbid.lower() + ".")), None)
            if file is not None:
                return file
            
def _get_cache_write_path(pdbid, compression = 'zst'):
    """Returns a file path where to write the downloaded PDB file."""
    if len(pdb_cache_dirs) == 0:
        return None
    pdb_cache_dir = pdb_cache_dirs[-1]
    subdir = os.path.join(pdb_cache_dir, pdbid[:2].lower())
    return os.path.join(subdir, pdbid.lower() + ".cif." + compression)

def _find_or_add_cache(pdbid):
    """Find a cif file in cache, or download it, add it to the cache and return the path."""
    cache_file = _find_in_cache(pdbid)
    if cache_file is not None:
        return cache_file
    else:
        write_path = _get_cache_write_path(pdbid, compression='zst')
        if write_path is None:
            raise Exception("No cache directory")
        os.makedirs(os.path.dirname(write_path), exist_ok=True)
        with download_pdb(pdbid) as f:
            write_path_ = write_path + ".part" + str(os.getpid())
            import zstandard
            with zstandard.open(write_path_, "wt", cctx=zstandard.ZstdCompressor(level=19)) as f2:
                f2.write(f.read())
            os.rename(write_path_, write_path)
        return write_path

def download_pdb(pdb_id: str) -> TextIO:
    """Downloads a CIF file from RCSB PDB"""
    import urllib.request
    # download from https://files.rcsb.org/download/XXXX.cif.gz
    url = f"https://files.rcsb.org/download/{pdb_id}.cif.gz"
    response = urllib.request.urlopen(url)
    return gzip.open(io.BytesIO(response.read()), "rt", encoding='utf-8-sig')

def get_pdb_file(file: Optional[str], pdbid: Optional[str] = None) -> str:
    if file is not None:
        return file
    return _find_or_add_cache(pdbid)
    
def open_pdb_file(file: Optional[str], pdbid: Optional[str] = None) -> TextIO:
    """
    Opens a PDBx CIF file as text reader.
    Either file or pdbid must be specified: open_pdb_file(None, '1ehz') or open_pdb_file('../data/1ehz.cif.gz')
    """
    if file is None and pdbid is not None:
        if len(pdb_cache_dirs) == 0:
            return download_pdb(pdbid)
        else:
            return open_pdb_file(_find_or_add_cache(pdbid), pdbid)

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
            print("Warning: unknown file extension", file)
            return open(file, "rt")
    elif file is None:
        raise Exception("No file specified")
    else:
        return file
    
def _gemmi_to_structure(cif):
    import gemmi
    if len(cif) != 1:
        raise Exception(f"CIF file has {len(cif)} blocks: {cif}")
    structure = gemmi.make_structure_from_block(cif[0])
    return structure, cif[0]

def load_pdb_gemmi(file: Optional[str | TextIO], pdbid: Optional[str] = None):
    """Loads a PDBx CIF file using gemmi library. Either file or pdbid must be specified: load_pdb_gemmi(None, '1ehz') or load_pdb_gemmi('../data/1ehz.cif.gz')"""
    import gemmi
    if isinstance(file, str):
        if file.endswith(".cif") or file.endswith(".cif.gz"):
            cif = gemmi.cif.read(file)
            return _gemmi_to_structure(cif)
    if isinstance(file, str) or file is None:
        with open_pdb_file(file, pdbid) as f:
            return load_pdb_gemmi(f, pdbid)
    else:
        data = file.read()
        cif = gemmi.cif.read_string(data)
        return _gemmi_to_structure(cif)
        # parser = Bio.PDB.MMCIFParser(QUIET=True)
        # structure = parser.get_structure(pdb_id, file)
        # return structure

def load_pdb(file: Optional[str | TextIO], pdbid: Optional[str] = None):
    """Loads a PDBx CIF file using BioPython library. Either file or pdbid must be specified: load_pdb(None, '1ehz') or load_pdb('../data/1ehz.cif.gz')"""
    import Bio.PDB.Structure
    if isinstance(file, str) or file is None:
        with open_pdb_file(file, pdbid) as f:
            return load_pdb(f, pdbid)
    else:
        try:
            parser = Bio.PDB.MMCIFParser(QUIET=True)
            structure = parser.get_structure(pdbid, file)
            return structure
        except Exception as e:
            print(f"Error loading {pdbid}: {e}")
            raise

@dataclasses.dataclass
class SymmetryOperation:
    """Symmetry operation specified in _pdbx_struct_oper_list PDB table"""
    # numeric? id
    id: str
    # 1_555 style code
    pdbname: str
    # x,y,z style code
    triplet: str
    # 3x3 rotation matrix
    rotation: np.ndarray
    translation: np.ndarray


@dataclasses.dataclass
class StructureData:
    """Some (meta)data extracted from mmCIF file."""
    assembly: list[SymmetryOperation]
    organism: Optional[str] = None
    organism_id: Optional[int] = None

def load_sym_data(file: Optional[str | TextIO], pdb_id: Optional[str] = None) -> StructureData:
    """
    Extracts some metadata from a mmCIF file. Either file or pdbid must be specified: load_sym_data(None, '1ehz') or load_sym_data('../data/1ehz.cif.gz')
    """
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

def load_sym_data_gemmi():
    import gemmi
    def core(b: gemmi.cif.Block):
        result = StructureData([])
        
        if obj := b.get_mmcif_category('_pdbx_struct_oper_list'):
            result.assembly = [
                SymmetryOperation(
                    id=obj['id'][i],
                    pdbname=obj['name'][i],
                    triplet=obj['symmetry_operation'][i],
                    rotation=np.array([
                        [ obj[f'matrix[1][1]'][i], obj[f'matrix[1][2]'][i], obj[f'matrix[1][3]'][i] ],
                        [ obj[f'matrix[2][1]'][i], obj[f'matrix[2][2]'][i], obj[f'matrix[2][3]'][i] ],
                        [ obj[f'matrix[3][1]'][i], obj[f'matrix[3][2]'][i], obj[f'matrix[3][3]'][i] ],
                    ]),
                    translation=np.array([ obj[f'vector[1]'][i], obj[f'vector[2]'][i], obj[f'vector[3]'][i] ])
                )
                for i in range(len(obj['id']))
            ]
        if obj := b.get_mmcif_category('_pdbx_entity_src_syn'):
            result.organism = obj['organism_scientific'][0]
            result.organism_id = int(obj['pdbx_src_id'][0]) if obj['pdbx_src_id'][0] else None
        return result

    return core
