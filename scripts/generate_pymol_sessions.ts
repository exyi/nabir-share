#!/usr/bin/env -S deno run --allow-read --allow-write --allow-run --allow-env
import { mergeReadableStreams } from "https://deno.land/std@0.195.0/streams/merge_readable_streams.ts";
import * as path from "https://deno.land/std/path/mod.ts";

const [ clustersDirectory, ] = Deno.args

const clusterFiles = Deno.statSync(clustersDirectory).isDirectory ? Array.from(Deno.readDirSync(clustersDirectory)).filter(d => d.isFile && d.name.endsWith(".csv")).map(d => clustersDirectory + "/" + d.name) : [ clustersDirectory ]

console.log(`Processing ${clusterFiles.join(", ")}`);

const script = path.resolve("./0_contacts_to_pml_both_frames")

const scriptProcesses = []

for (const clusterFile of clusterFiles) {
    const cwd = path.dirname(clusterFile)
    const clusterName = path.basename(clusterFile, ".csv")
    const pymolScript = new Deno.Command(script, { cwd, args: [ clusterName + ".csv" ], stdout: "piped", stdin: "null" }).spawn()
    pymolScript.stdout.pipeTo((await Deno.open(cwd + "/" + clusterName + ".pml", { create: true, write: true })).writable)

    scriptProcesses.push(pymolScript)
}

await Promise.all(scriptProcesses.map(p => p.status))
console.log(`Created ${scriptProcesses.length} pymol scripts, running pymol on all of them...`)

const pymolProcesses = []
for (const clusterFile of clusterFiles) {
    const cwd = path.dirname(clusterFile)
    const clusterName = path.basename(clusterFile, ".csv")
    const pymolProcess = new Deno.Command("pymol", { cwd, args: [ "-c", clusterName + ".pml" ], stdout: "piped", stderr: "piped", stdin: "null" }).spawn()

    const outputStream = mergeReadableStreams(pymolProcess.stdout, pymolProcess.stderr)
    outputStream.pipeTo((await Deno.open(cwd + "/" + clusterName + ".pml.log", { create: true, write: true })).writable)
    pymolProcesses.push({ pymolProcess, clusterFile })
}

await Promise.all(pymolProcesses.map(async p => {
    await p.pymolProcess.status
    console.log(`Finished running pymol on ${p.clusterFile}`)
}))
