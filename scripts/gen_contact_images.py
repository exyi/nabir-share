import pymol
from pymol import cmd
import polars as pl
import os, sys
import pdb_utils


def residue_selection(chain, nt, ins, alt):
    chain = str(chain)
    nt = str(nt).replace("-", "\\-")
    if ins and ins != ' ' and ins != '?':
        nt += str(ins)

    if alt and alt != ' ' and alt != '?':
        alt = f" alt {alt}"
    else:
        alt = ""
    
    return f"(chain {chain} and resi {nt}{alt})"

def load(file, pdbid):
    if file:
        cmd.load(file)
    else:
        if len(pdb_utils.pdb_cache_dirs) == 0:
            cmd.fetch(pdbid)
        else:
            file = pdb_utils.get_pdb_file(None, pdbid)
            print(file)
            cmd.load(file)

def orient_pair(pdbid, chain1, nt1, ins1, alt1, chain2, nt2, ins2, alt2):
    cmd.hide("everything", f"%{pdbid}")
    pair_selection =f"%{pdbid} and ({residue_selection(chain1, nt1, ins1, alt1)} or {residue_selection(chain2, nt2, ins2, alt2)})"
    print(pair_selection)
    cmd.select("pair", pair_selection)
    cmd.show("sticks", "%pair")
    cmd.distance("pair_contacts", "%pair", "%pair", mode=2)
    cmd.orient("%pair")
    cmd.util.cbag(f"%pair")
    cmd.show("lines", f"%{pdbid}")
    cmd.label("%pair", "name")
    cmd.clip("slab", 12, "%pair")
    cmd.color("grey", f"%{pdbid} and not %pair")

def make_pair_image(output_file, pdbid, chain1, nt1, ins1, alt1, chain2, nt2, ins2, alt2):
    orient_pair(pdbid, chain1, nt1, ins1, alt1, chain2, nt2, ins2, alt2)
    cmd.png(output_file, width=640, height=640, ray=1)
    print(f"Saved {output_file}")

def process_group(pdbid, group: pl.DataFrame, output_dir: str):
    pdbid = str(pdbid)
    load(None, pdbid)
    pdbdir = os.path.join(output_dir, pdbid)
    os.makedirs(pdbdir, exist_ok=True)
    for chain1, nt1, ins1, alt1, chain2, nt2, ins2, alt2 in zip(group["chain1"], group["nr1"], group["ins1"], group["alt1"], group["chain2"], group["nr2"], group["ins2"], group["alt2"]):
        ins1 = None if ins1 == ' ' else ins1
        ins2 = None if ins2 == ' ' else ins2
        alt1 = None if alt1 == '?' else alt1
        alt2 = None if alt2 == '?' else alt2
        output_file = os.path.join(pdbdir, f"{chain1}_{nt1}{ins1 or ''}{alt1 or ''}-{chain2}_{nt2}{ins2 or ''}{alt2 or ''}.png")
        make_pair_image(output_file, pdbid, chain1, nt1, ins1, alt1, chain2, nt2, ins2, alt2)

def make_images(df: pl.DataFrame, output_dir: str):
    for pdbid, group in sorted(df.groupby("pdbid")):
        process_group(pdbid, group, output_dir)

def make_images_mp(df: pl.DataFrame, output_dir: str, threads: int):
    import multiprocessing
    multiprocessing.set_start_method("spawn")
    with multiprocessing.Pool(threads) as pool:
        processes = [
            pool.apply_async(process_group, (pdbid, group, output_dir))

            for pdbid, group in sorted(df.groupby("pdbid"))
        ]
        for p in processes:
            p.get()

def main(argv):
    import argparse
    parser = argparse.ArgumentParser(description="Generate contact images")
    parser.add_argument("input", help="Input CSV file", nargs="+")
    parser.add_argument("--output-dir", "-o", required=True, help="Output directory")
    parser.add_argument("--threads", "-t", type=int, default=1, help="Number of threads to use")
    args = parser.parse_args(argv)

    import pair_csv_parse
    df = pair_csv_parse.scan_pair_csvs(args.input).collect()
    if args.threads == 1:
        make_images(df, args.output_dir)
    else:
        make_images_mp(df, args.output_dir, args.threads)

if __name__ == "__main__" or __name__ == "pymol":
    main(sys.argv[1:])


