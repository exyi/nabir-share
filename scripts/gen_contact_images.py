import itertools
import tempfile
from typing import Optional
import pymol
from pymol import cmd
import polars as pl
import os, sys
import pdb_utils
import numpy as np
import math
from dataclasses import dataclass
import shutil
import subprocess
import pair_defs

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

def transform_to_camera_space(coords):
    coords = np.array(coords)
    view = cmd.get_view()
    matrix = np.array(view[0:9]).reshape((3, 3))
    camera_center = np.array(view[9:12])
    model_center = np.array(view[12:15])

    camera_coords = camera_center + np.dot(matrix.T, coords - model_center)
    return camera_coords

def rotate_to_y_axis(bond1, bond2):
    coord1 = transform_to_camera_space(cmd.get_coords(bond1)[0])
    coord2 = transform_to_camera_space(cmd.get_coords(bond2)[0])
    # rotate around z-axis such that coord1 is right above coord2
    angle = np.arctan2(coord2[0] - coord1[0], coord2[1] - coord1[1])
    angle = angle / math.pi * 180
    # print("the bond coordinate is ", coord1, coord2)
    # print(f"1: {np.arctan2(coord1[1], coord1[0]) / math.pi * 180}, 2: {np.arctan2(coord2[1], coord2[0]) / math.pi * 180}")
    print("rotating by ", angle)
    cmd.turn("z", angle)
    # cmd.color("red", f"({bond1}) or ({bond2})")
    return abs(angle) > 0.1

def orient_nucleotide_as_main():
    if cmd.get_coords("%rightnt and (name N9)") is None:
        natom = "N1"
    else:
        natom = "N9"
    for i in range(1):
        # it does not converge instantly, no idea why
        rotate_to_y_axis("%rightnt and (name C1')", f"%rightnt and (name {natom})")

    # main nucleotide should be on the left
    rightC1_A = transform_to_camera_space(cmd.get_coords("%rightnt and (name C1')")[0])
    leftC1_A = transform_to_camera_space(cmd.get_coords("%pair and (not %rightnt) and (name C1')")[0])
    if rightC1_A[0] > leftC1_A[0]:
        cmd.turn("y", 180)
        rightC1_B = transform_to_camera_space(cmd.get_coords("%rightnt and (name C1')")[0])
        print(f"flipped {rightC1_A} to {rightC1_B}")
        rotate_to_y_axis("%rightnt and (name C1')", f"%rightnt and (name {natom})")

def find_and_select_water(nt1, nt2):
    # coords1 = np.array(cmd.get_coords(nt1))
    # coords2 = np.array(cmd.get_coords(nt2))
    # coords_water = np.array(cmd.get_coords("resn HOH and name O"))
    # dist1 = np.min(np.linalg.norm(coords_water.reshape(1, -1, 3) - coords1.reshape(-1, 1, 3), axis=2), axis=0)
    # assert dist1.shape == (coords_water.shape[0],)
    # dist2 = np.min(np.linalg.norm(coords_water.reshape(1, -1, 3) - coords2.reshape(-1, 1, 3), axis=2), axis=0)
    # assert dist2.shape == (coords_water.shape[0],)

    threshold = 3.6

    cmd.select("nwaters", f"(resn HOH within {threshold} of ({nt1})) and (resn HOH within {threshold} of ({nt2}))")

def get_margins(image):
    cmd.png(image, width=80, height=45, ray=1)#, quiet=1)
    import imageio.v3 as iio
    img = iio.imread(image)
    img = np.any(img[:, :, 0:3] > 0, axis=2)
    imgx = np.any(img, axis=0)
    imgy = np.any(img, axis=1)
    if not np.any(imgx):
        return None

    return np.array([
        np.argmax(imgx),
        np.argmax(imgy),
        np.argmax(imgx[::-1]),
        np.argmax(imgy[::-1])
    ])

def zoom_to_borders(center):
    with tempfile.NamedTemporaryFile(suffix=".png") as f:
        observations = []
        def experiment(zoom):
            cmd.zoom(center, zoom)
            m = get_margins(f.name)
            observations.append((zoom, np.min(m) if m is not None else 0, m))
            return m

        m = experiment(0)
        if m is None:
            print("WARNING: empty image!")
            return
        while True:
            best = min((o for o in observations if o[1] > 0), key=lambda x: x[1])

def highlight_bonds(res1, res2, hbonds: Optional[list[tuple[str, str, str, str]]], label_atoms):
    if hbonds is None:
        cmd.distance("pair_contacts", res1, res2, mode=2)
        return
    
    def atom_sele(atom: str):
        if atom.startswith("A"):
            return f"({res1}) and name \"{atom[1:]}\""
        elif atom.startswith("B"):
            return f"({res2}) and name \"{atom[1:]}\""
        else:
            return f"name \"{atom}\""

    labels = []
    for i, (atom0, atom1, atom2, atom3) in enumerate(hbonds):
        cmd.distance(f"pair_hbond_{i}", atom_sele(atom1), atom_sele(atom2), mode=0)
        cmd.hide("labels", f"pair_hbond_{i}")
        if label_atoms:
            labels.extend([atom1, atom2])

    if len(labels) > 0:
        cmd.label("(" + " or ".join([
            "(" + atom_sele(a) + ")"
            for a in labels
        ]) + ")", "name")

def orient_pair(pdbid, chain1, nt1, ins1, alt1, chain2, nt2, ins2, alt2,
    standard_orientation,
    hbonds: Optional[list[tuple[str, str, str, str]]] = None,
    label_atoms=True,
    grey_context=False,
    find_water=False):
    cmd.hide("everything", f"%{pdbid}")
    pair_selection =f"%{pdbid} and ({residue_selection(chain1, nt1, ins1, alt1)} or {residue_selection(chain2, nt2, ins2, alt2)})"
    print(pair_selection)
    cmd.select("pair", pair_selection)
    cmd.select("rightnt", f"%pair and ({residue_selection(chain1, nt1, ins1, alt1)})")
    # cmd.select("rightnt", f"%pair and ({residue_selection(chain2, nt2, ins2, alt2)})")
    cmd.show("sticks", "%pair")
    cmd.delete("pair_contacts")
    cmd.delete("pair_w_contacts")
    for i in range(5):
        cmd.delete(f"pair_hbond_{i}")
    # cmd.distance("pair_contacts", "%pair", "%pair", mode=2)
    # cmd.hide("labels", "pair_contacts")
    if find_water:
        find_and_select_water("%rightnt", "%pair and not %rightnt")
        cmd.distance("pair_w_contacts", "%pair", "%nwaters", mode=2)
        cmd.hide("labels", "pair_w_contacts")
        cmd.show("nb_spheres", "%nwaters")
        cmd.color("0xff0000", "%nwaters")
        # cmd.set("nb_spheres_size", 0.3, "nwaters")
    cmd.util.cba("grey", f"%pair")
    normal_atoms_selection = "not (name C2' or name C3' or name C4' or name C5' or name O2' or name O3' or name O4' or name O5' or name P or name OP1 or name OP2)"
    highlight_bonds("%rightnt", "%pair and not %rightnt", hbonds, label_atoms)
    # if "all" in label_atoms:
    #     cmd.label("%pair", "name")
    if standard_orientation:
        cmd.orient(f"%rightnt and {normal_atoms_selection}")
        orient_nucleotide_as_main()
    else:
        cmd.orient(f"%pair and {normal_atoms_selection}")
    cmd.zoom("%pair", 0)
    # cmd.center("%rightnt")
    cmd.set("label_color", "black")
    cmd.set("label_bg_color", "white", "%pair")
    cmd.set("label_outline_color", "white")
    cmd.set("label_size", 25)
    # cmd.set("label_position", (0, 0, 0))
    cmd.set("label_font_id", 7)
    if grey_context:
        cmd.show("sticks", f"%{pdbid} and not %pair")
        cmd.set_bond("stick_radius", 0.07, f"%{pdbid} and not %pair and not resname HOH")
        cmd.set_bond("stick_transparency", 0.4, f"%{pdbid} and not %pair and not resname HOH")
        cmd.clip("slab", 12, "%pair")
        cmd.color("0xeeeeee", f"%{pdbid} and not %pair and not resname HOH")

    cmd.set_bond("stick_radius", 0.25, "%pair")
    cmd.set_bond("stick_transparency", 0.0, "%pair")

@dataclass
class BPArgs:
    # label_atoms: list[str]
    standard_orientation: bool
    movie: int
    incremental: bool
    skip_bad: bool

def make_pair_image(output_file, pdbid, chain1, nt1, ins1, alt1, chain2, nt2, ins2, alt2, bpargs: BPArgs, hbonds: Optional[list[tuple[str, str, str, str]]]):
    print(bpargs)
    orient_pair(pdbid, chain1, nt1, ins1, alt1, chain2, nt2, ins2, alt2, standard_orientation=bpargs.standard_orientation, hbonds=hbonds)
    cmd.png(output_file, width=2560, height=1440, ray=1)
    print(f"Saved basepair image {output_file}")

def make_pair_rot_movie(output_file, pdbid, chain1, nt1, ins1, alt1, chain2, nt2, ins2, alt2, bpargs: BPArgs, hbonds: Optional[list[tuple[str, str, str, str]]]):
    length = bpargs.movie
    if length == 0:
        return

    orient_pair(pdbid, chain1, nt1, ins1, alt1, chain2, nt2, ins2, alt2, standard_orientation=bpargs.standard_orientation, hbonds=hbonds, label_atoms=False, grey_context=True, find_water=True)
    try:
        cmd.mset("1", length)
        cmd.mview("store", 1)
        cmd.mview("store", length)
        cmd.turn("x", 120)
        cmd.mview("store", length // 3, power=1)
        cmd.turn("x", 120)
        cmd.mview("store", length - length // 3, power=1)
        cmd.bg_color("white")
        cmd.set("ray_shadow", 0)
        cmd.set("antialias", 1)
        # cmd.set("max_threads", 1!)

        use_ffmpeg_directly = True
        if output_file.endswith(".dir"):
            os.makedirs(output_file, exist_ok=True)
            cmd.mpng(output_file + "/", width=640, height=480)
        elif use_ffmpeg_directly:
            os.makedirs(output_file + ".pngdir", exist_ok=True)
            try:
                # cmd.mpng(output_file + ".pngdir/", mode=2, width=1920, height=1080)
                cmd.mpng(output_file + ".pngdir/", mode=2, width=1280, height=720)
                if os.path.exists(output_file):
                    os.remove(output_file)
                ffmpeg_result = subprocess.run([
                    "ffmpeg",
                    "-framerate", "30",
                    "-pattern_type", "glob",
                    "-i", output_file + ".pngdir/*.png",
                    "-c:v", "libx264",
                    "-pix_fmt", "yuv420p",
                    output_file
                ], capture_output=True)
                if ffmpeg_result.returncode != 0:
                    print("ffmpeg failed")
                    print(ffmpeg_result.stdout.decode("utf-8"))
                    print(ffmpeg_result.stderr.decode("utf-8"))
            finally:
                shutil.rmtree(output_file + ".pngdir")
        else:
            from pymol import movie
            movie.produce(output_file, mode="draw", quality=96, width=1280, height=720, encoder="ffmpeg")
        print("movie produced")
    finally:
        cmd.mclear()


def process_group(pdbid, group: pl.DataFrame, output_dir: str, bpargs: BPArgs):
    if bpargs.skip_bad:
        orig_len = len(group)
        bond_length_cols = [pl.col(f"hb_{i}_length") for i in range(5) if f"hb_{i}_length" in group.columns]
        if len(bond_length_cols):
            group = group.filter(
                pl.any_horizontal(*[c.is_not_null() & (c < 5) for c in bond_length_cols])
            )
        print(f"skipping {orig_len - len(group)}/{orig_len} bad pairs in {pdbid}")

    if len(group) == 0:
        return
    pdbid = str(pdbid)
    loaded = False
    def load_if_needed() -> None:
        nonlocal loaded
        if not loaded:
            load(None, pdbid)
            loaded = True
    pdbdir = os.path.join(output_dir, pdbid)
    os.makedirs(pdbdir, exist_ok=True)
    type_column = group["type"] if "type" in group.columns else itertools.repeat(None)
    for pair_type, res1, chain1, nt1, ins1, alt1, res2, chain2, nt2, ins2, alt2 in zip(type_column, group["res1"], group["chain1"], group["nr1"], group["ins1"], group["alt1"], group["res2"], group["chain2"], group["nr2"], group["ins2"], group["alt2"]):
        ins1 = None if ins1 == ' ' else ins1
        ins2 = None if ins2 == ' ' else ins2
        alt1 = None if alt1 == '?' else alt1
        alt2 = None if alt2 == '?' else alt2
        if pair_type:
            pt = pair_defs.PairType(pair_type, (pair_defs.map_resname(res1), pair_defs.map_resname(res2)))
            hbonds = pair_defs.get_hbonds(pt)
            hbonds = [ b for b in hbonds if not pair_defs.is_bond_hidden(pt, b) ]
        else:
            hbonds = None

        output_file = os.path.join(pdbdir, f"{chain1}_{nt1}{ins1 or ''}{alt1 or ''}-{chain2}_{nt2}{ins2 or ''}{alt2 or ''}")
        if bpargs.incremental and os.path.exists(output_file + ".png"):
            if bpargs.movie == 0 or os.path.exists(output_file + ".mp4"):
                print(f"Skipping {output_file} because it already exists")
                continue

        load_if_needed()
        make_pair_image(output_file + ".png", pdbid, chain1, nt1, ins1, alt1, chain2, nt2, ins2, alt2, bpargs, hbonds)
        make_pair_rot_movie(output_file + ".mp4", pdbid, chain1, nt1, ins1, alt1, chain2, nt2, ins2, alt2, bpargs, hbonds)

    cmd.reinitialize()

def make_images(df: pl.DataFrame, output_dir: str, bpargs: BPArgs):
    for pdbid, group in sorted(df.groupby("pdbid")):
        process_group(pdbid, group, output_dir, bpargs)
    # if cmd.is_gui_thread():
        # cmd.quit()

def process_init(niceness, affinity):
    import multiprocessing
    if niceness:
        x = os.nice(0)
        os.nice(niceness - x)
    if affinity:
        os.sched_setaffinity(0, affinity)

def make_images_mp(df: pl.DataFrame, output_dir: str, threads: int, affinity: Optional[list[int]], niceness: Optional[int], bpargs: BPArgs):
    # if affinity is None:
    #     affinities = [None] * threads
    # elif len(affinity) == threads:
    #     affinities = [ [a] for a in affinity ]
    # else:
    #     affinities = [affinity] * threads
    import multiprocessing
    multiprocessing.set_start_method("spawn")
    with multiprocessing.Pool(threads, initializer=process_init, initargs=(niceness,affinity)) as pool:
        processes = [
            pool.apply_async(process_group, (pdbid, group, output_dir, bpargs))

            for pdbid, group in sorted(df.groupby("pdbid"))
        ]
        for p in processes:
            p.get()

def main(argv):
    import argparse
    parser = argparse.ArgumentParser(description="Generate contact images")
    parser.add_argument("input", help="Input CSV file", nargs="+")
    parser.add_argument("--output-dir", "-o", required=True, help="Output directory")
    parser.add_argument("--threads", "-t", type=int, default=None, help="Number of threads to use")
    parser.add_argument("--cpu_affinity", type=int, nargs="*", default=None, help="Which CPUs to use (list of integers, indexed from 0)")
    parser.add_argument("--niceness", type=int, default=None, help="Run the process with the specified niceness (don't change by default)")
    parser.add_argument("--standard-orientation", type=eval, default=False, help="When set to true, orient the base pair such that the first nucleotide is always left and the N1/N9 - C1' is along the y-axis (N is above C)")
    parser.add_argument("--label-atoms", type=str, nargs="*", default=[], help="Atom names to label")
    parser.add_argument("--movie", type=int, default=0, help="If not zero, produce a rotating animation of the base pair")
    parser.add_argument("--incremental", type=bool, default=False, help="Generate the image/video only if it does not exist yet")
    parser.add_argument("--skip-bad", type=bool, default=False, help="Skip basepairs with bad or missing parameters")
    args = parser.parse_args(argv)

    import pair_csv_parse
    df = pair_csv_parse.scan_pair_csvs(args.input).collect()
    bpargs = BPArgs(args.standard_orientation, args.movie, args.incremental, args.skip_bad)
    threads = args.threads
    if threads is None and args.cpu_affinity is not None:
        threads = len(args.cpu_affinity)
    if threads is None:
        threads = 1

    if threads == 1:
        process_init(args.niceness, args.cpu_affinity)
        make_images(df, args.output_dir, bpargs)
    else:
        make_images_mp(df, args.output_dir, args.threads, args.cpu_affinity, args.niceness, bpargs)

if __name__ == "__main__" or __name__ == "pymol":
    main(sys.argv[1:])


