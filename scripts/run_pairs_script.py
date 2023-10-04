#!/usr/bin/env python3

import os, sys, re, argparse, subprocess, threading, time, math

def run(outfile, input_files, pdbcache, threads, pair_type):
    a = [
        "python3", "./pairs.py",
        *input_files,
        "--pdbcache", *pdbcache,
        "--threads", str(threads),
        "--pair-type", pair_type,
        "--dedupe",
        "--output", outfile + ".csv",
    ]
    print("Running pairs: " + outfile)
    print("Running command: " + " ".join(a))
    starttime = time.time()
    proc = subprocess.Popen(a, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    outlock = threading.Lock()
    with open(outfile + ".log", "wt") as f:
        def write(stream, streamname):
            for line in stream:
                with outlock:
                    t = time.time() - starttime
                    f.write(f"[{t:9.2f} {streamname}] " + line.decode("utf-8"))
        t_out = threading.Thread(target=write, args=[proc.stdout, "OUT"])
        t_out.start()
        t_err = threading.Thread(target=write, args=[proc.stderr, "ERR"])
        t_err.start()
        t_out.join()
        t_err.join()
    
    if proc.wait() != 0:
        print("Error running pairs: " + outfile)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run pairs script")
    parser.add_argument("input_directory", help="Directory with cWW/A_G.csv files")
    parser.add_argument("pairs", help="List of pair types to analyze, format: G-C-cWW", nargs="*")
    parser.add_argument("--pdbcache", nargs="+", help="Directories to search for PDB files in order to avoid downloading. Last directory will be written to, if the structure is not found and has to be downloaded.")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use")
    parser.add_argument("--output-dir", "-o", required=True, help="Output directory")
    args = parser.parse_args()

    ps = args.pairs
    if ps == ["all"]:
        import pairs
        ps = [ f"{b}-{t}" for t, b in pairs.hbonding_atoms.keys() ]
    our_threads = max(1, min(len(ps) // 2, args.threads // 2, 4))
    sub_threads = round(args.threads / our_threads)
    semaphore = threading.Semaphore(our_threads)

    def thread_main(pair: str):
        semaphore.acquire()
        try:
            pair_base1, pair_base2, pair_type = pair.split("-")
            input_dir = os.path.join(args.input_directory, pair_type)
            file_regex = f"D?{re.sub('[TU]', '[TU]', pair_base1)}_D?{re.sub('[TU]', '[TU]', pair_base2)}[.]csv"
            in_files = [ os.path.join(input_dir, f) for f in os.listdir(input_dir) if re.match(file_regex, f) ]
            in_files = [ f for f in in_files if os.path.isfile(f) and os.stat(f).st_size > 0 ]
            assert len(in_files) > 0, f"Could not find any files matching {file_regex} in {input_dir}"
            run(outfile=os.path.join(args.output_dir, pair),
                input_files=in_files,
                pdbcache=args.pdbcache,
                threads=sub_threads,
                pair_type=pair_type
            )
        finally:
            semaphore.release()

    threads = [ threading.Thread(target=thread_main, args=[pair]) for pair in ps ]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
