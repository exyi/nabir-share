<script lang="ts">
  import { initDB } from '$lib/duckdb'
  import PairImages from '$lib/components/pairimages.svelte'
  import Spinner from '$lib/components/Spinner.svelte'
  export let data;
  import { base } from '$app/paths';
  import metadata from '$lib/metadata'
	import FilterEditor from '$lib/components/filterEditor.svelte';
	import { filterToSqlCondition, makeSqlQuery, type NucleotideFilterModel } from '$lib/dbLayer.js';
	import { fix_position } from 'svelte/internal';
	import { parsePairingType, type NucleotideId, type PairId, type PairingInfo } from '$lib/pairing.js';
  const parquetFiles = {
    'A-G-tSS': `/pairing_tables/A-G-tSS.csv.parquet`,
    'A-G-tHS': `/pairing_tables/A-G-tHS.csv.parquet`,
    'G-G-cWH': `/pairing_tables/G-G-cWH.csv.parquet`,
  }
  async function load_db() {
    console.log("LOADING DB")
    // A simple case of a db of a single parquet file.
    const db = await initDB();
    for (const [name, url] of Object.entries(parquetFiles)) {
      await db.registerFileURL(name, `${base}${url}`, 4, false);
    }
    const conn = await db.connect();

    console.log("PREPARING VIEWS")
    for (const [name, url] of Object.entries(parquetFiles)) {
      await conn.query(`CREATE VIEW '${name}' AS SELECT * FROM parquet_scan('${name}')`);
    }
    // await conn.query(`CREATE TABLE p1 AS SELECT * FROM parquet_scan('SOTU.parquet')`);
    // await conn.query(`CREATE VIEW wordcounts_raw AS SELECT * FROM (SELECT "@id" id, 
    //     UNNEST("nc:unigrams").word0 as word, 
    //     UNNEST("nc:unigrams").count as count FROM
    //     p1) t1
    // `);
    // await conn.query(`
    //   CREATE TABLE wordcounts AS SELECT * FROM 
    //   wordcounts_raw
    //   NATURAL JOIN (SELECT word, SUM(count) as tot, COUNT(*) AS df FROM wordcounts_raw GROUP BY word) t2
    // `);
    return conn;
  }

  let selectedPairing = 'G-G-cWH'
  let filter: NucleotideFilterModel = { bond_acceptor_angle: [], bond_donor_angle: [], bond_length: [] }

  // Set up the db connection as an empty promise.
  const conn_prom = load_db();

  let resultsPromise = new Promise(() => {})
  let results = []

  $: {
    filter, selectedPairing
    updateResults()
  }

  function* convertQueryResults(rs, pairingType, limit=undefined): Generator<PairingInfo> {
    function convName(name) {
      if (name == null) return undefined
      name = String(name).trim()
      if (name == '' || name == '?' || name == '\0') return undefined
      return name
    }

    let c = 0
    for (const r of rs) {
      if (limit != null && c >= limit)
        break
      const pdbid = r.pdbid, model = Number(r.model)
      const nt1: NucleotideId = { pdbid, model, chain: convName(r.chain1), resnum: Number(r.nr1), resname: convName(r.res1), altloc: convName(r.alt1), inscode: convName(r.ins1) }
      const nt2: NucleotideId = { pdbid, model, chain: convName(r.chain2), resnum: Number(r.nr2), resname: convName(r.res2), altloc: convName(r.alt2), inscode: convName(r.ins2) }

      const id: PairId = { nt1, nt2, pairingType }
      yield { id }
      c++
    }
  }

  async function updateResults() {
    const conn = await conn_prom;
    const sql = makeSqlQuery(filter, `parquet_scan('${selectedPairing}')`)
    console.log(sql)
    resultsPromise = conn.query(sql)
    results = Array.from(convertQueryResults(await resultsPromise as any, parsePairingType(selectedPairing), 100));
    console.log({resultsPromise, results})
  }

  // const imgDir = "http://[2a01:4f8:c2c:bb6c:4::8]:12345/"
  const imgDir = "https://pairs.exyi.cz/img"
  // const imgDir = base+"/img"
</script>


<h1>Nucleotide base pairing visualizer</h1>

<div class="selector">
  {#each Object.keys(parquetFiles) as p}
    <button class:selected={selectedPairing == p} on:click={() => {selectedPairing = p}}>{p}</button>
  {/each}
</div>
<div class="filters">
  <FilterEditor bind:filter={filter} isPlainSql={false} />
</div>
{#await resultsPromise}
<div style="position: relative">
  <div style="position: absolute; top: 0; left: 0">
    <Spinner></Spinner>
  </div>
</div>
{:then result}
{/await}
<PairImages pairs={results} rootImages={imgDir} imgAttachement=".png" videoAttachement=".mp4" />
<style>
.selector .selected {
  border: 2px solid black;
}
.selector {

}
</style>
