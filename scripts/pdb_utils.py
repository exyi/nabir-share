
import os, gzip, io
from typing import Any, Callable, Optional, TextIO, TypeVar, Union


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
        subdir = os.path.join(cache_dir, pdb_id[:2].lower())
        if not os.path.isdir(subdir):
            subdir = ([ os.path.join(cache_dir, s) for s in os.listdir(cache_dir) if pdb_id.lower().startswith(s.lower())] or [ None ])[0]

        if subdir is not None:
            file = ([ os.path.join(subdir, s) for s in os.listdir(subdir) if s.lower().startswith(pdb_id.lower() + ".") ] or [None])[0]
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

