<script lang="ts">
  import PairImages from '$lib/components/pairimages.svelte'
  import Spinner from '$lib/components/Spinner.svelte'
  export let data;
  import { base } from '$app/paths';
  import metadata from '$lib/metadata'
	import FilterEditor from '$lib/components/filterEditor.svelte';
	import { aggregateBondParameters, aggregatePdbCountQuery, aggregateTypesQuery, filterToSqlCondition, makeSqlQuery, type NucleotideFilterModel } from '$lib/dbLayer.js';
	import { fix_position } from 'svelte/internal';
	import { parsePairingType, type NucleotideId, type PairId, type PairingInfo, type HydrogenBondInfo, tryParsePairingType } from '$lib/pairing.js';
	import { Modal } from 'svelte-simple-modal';
	import type { AsyncDuckDBConnection } from '@duckdb/duckdb-wasm';
	import { AsyncLock } from '$lib/lock.js';
  import * as db from '$lib/dbInstance';
  import type * as arrow from 'apache-arrow'

  let selectedFamily = 'tWW'
  let selectedPairing = 'tWW-A-A'

  let filterMode: "ranges" | "sql" = "ranges"
  let filter: NucleotideFilterModel = { bond_acceptor_angle: [], bond_donor_angle: [], bond_length: [], filtered: true, includeNears: false }

  // Set up the db connection as an empty promise.
  const connPromise = db.connect();

  let tableLoadPromise
  $: {
    selectedPairing;

    onSelectedPairingChange()
  }


  let lastSelectedPairing = null
  function onSelectedPairingChange() {
    if (lastSelectedPairing == selectedPairing)
      return
    lastSelectedPairing = selectedPairing

    tableLoadPromise = updateResultsLock.withLock(() =>
      connPromise
        .then(async conn => {
          console.log("Dropping existing views, switching to ", selectedPairing)
          // await conn.cancelSent()
          await conn.query(`DROP VIEW IF EXISTS selectedpair`)
          await conn.query(`DROP VIEW IF EXISTS selectedpair_f`)
          await conn.query(`DROP VIEW IF EXISTS selectedpair_n`)
          // return Promise.all([
          //   conn.query(`CREATE OR REPLACE VIEW 'selectedpair' AS SELECT * FROM parquet_scan('${selectedPairing}')`),
          //   conn.query(`CREATE OR REPLACE VIEW 'selectedpair_f' AS SELECT * FROM parquet_scan('${selectedPairing}-filtered')`),
          //   // conn.query(`CREATE OR REPLACE VIEW 'selectedpair_n' AS SELECT * FROM parquet_scan('${selectedPairing}-n')`)
          // ])
        })
    )

    tableLoadPromise.catch(e => {
        console.error(`Could not load table ${selectedPairing}:`, e)
    })
  }

  async function ensureViews(conn: AsyncDuckDBConnection, abort: AbortSignal, query: string) {
    async function addView(name: string, file: string) {
      console.log(`Loading ${file} into view ${name}`)
      abort.throwIfAborted()
      await conn.query(`CREATE OR REPLACE VIEW '${name}' AS SELECT * FROM parquet_scan('${file}')`)
    }
    // await conn.
    const queryTables = new Set(await conn.getTableNames(query))
    const existingTables = new Set([...await conn.query("select view_name as name from duckdb_views() union select table_name as name from duckdb_tables()")].map(r => r.name))

    for (const e in existingTables) {
      queryTables.delete(e)
    }
    console.log("missing tables:", [...queryTables])
    if (queryTables.has('selectedpair')) {
      await addView('selectedpair', `${selectedPairing}`)
      queryTables.delete('selectedpair')
    }
    if (queryTables.has('selectedpair_f')) {
      await addView('selectedpair_f', `${selectedPairing}-filtered`)
      queryTables.delete('selectedpair_f')
    }
    if (queryTables.has('selectedpair_n')) {
      await addView('selectedpair_n', `n${selectedPairing}`)
      queryTables.delete('selectedpair_n')
    }

    for (const t of queryTables) {
      if (tryParsePairingType(t) != null) {
        await addView(t, t)
      } else if ((t.endsWith('_f') || t.endsWith('-f') || t.endsWith('_n') || t.endsWith('-n')) &&
        tryParsePairingType(t.slice(0, -2)) != null) {
        await addView(t, (t.endsWith('n') ? 'n' : '') + t.slice(0, -2) + (t.endsWith('f') ? '-filtered' : ''))
      } else {
        if (!existingTables.has(t))
          console.warn(`Maybe missing table: ${t}?`)
      }
    }
  }

  type ResultsAggregates = {
    types?: { [key: string]: number },
    pdbStructures?: { [key: string]: number },
    bondStats?: { bond: number, stat: "nncount" | "mean" | "min" | "max" | "median" | "p10" | "p25" | "p75" | "p90" | "stddev", param: "length" | "acceptor_angle" | "donor_angle", value: number }[]
  }
  let resultsPromise: Promise<arrow.Table<any> | undefined> = new Promise(() => {})
  let results = []
  let resultsCount = 0
  let resultsAgg: ResultsAggregates = {}

  $: {
    filter, filterMode, selectedPairing
    updateResults()
  }

  const requiredColumns = [ "pdbid", "chain1", "nr1", "chain2", "nr2", ]
  const recommendedColumns = [ "model", "ins1", "alt1", "res1", "res2", "ins2", "alt2", "res2" ]

  function getMetadata(pairType) {
    return metadata.find(m => m.pair_type.join('-').toLowerCase() == pairType.toLowerCase())
  }

  function* convertQueryResults(rs: arrow.Table<any>, pairingType, limit=undefined): Generator<PairingInfo> {
    function convName(name) {
      if (name == null) return undefined
      name = String(name).trim()
      if (name == '' || name == '?' || name == '\0') return undefined
      return name
    }

    const meta = getMetadata(selectedPairing)

    let c = 0
    for (const r of rs) {
      if (limit != null && c >= limit)
        break
      const pdbid = r.pdbid, model = Number(r.model ?? 1)
      const nt1: NucleotideId = { pdbid, model, chain: convName(r.chain1), resnum: Number(r.nr1), resname: convName(r.res1), altloc: convName(r.alt1), inscode: convName(r.ins1) }
      const nt2: NucleotideId = { pdbid, model, chain: convName(r.chain2), resnum: Number(r.nr2), resname: convName(r.res2), altloc: convName(r.alt2), inscode: convName(r.ins2) }

      const id: PairId = { nt1, nt2, pairingType }
      let hbonds: HydrogenBondInfo[] | undefined
      if ("hb_0_length" in r) {
        hbonds = []
        for (let i = 0; i <= 4; i++) {
          const length = r[`hb_${i}_length`]
          if (length == null) break
          hbonds.push({
            length: Number(length),
            donorAngle: Number(r[`hb_${i}_donor_angle`]),
            acceptorAngle: Number(r[`hb_${i}_acceptor_angle`]),
            label: meta?.labels[i],
          })
        }
      }
      const coplanarity = r.bogopropeller ?? r.coplanarity
      yield { id, hbonds, coplanarity, originalRow: r }
      c++
    }
  }

  const updateResultsLock = new AsyncLock()

  function queryFromTable(filter: NucleotideFilterModel) {
    let queryTable = `selectedpair` + (filter.filtered ? "_f" : "")
    if (filter.includeNears) {
      queryTable = `(select * FROM ${queryTable} UNION ALL BY NAME SELECT * from selectedpair_n)`
    }
    return queryTable
  }
  async function updateResults() {
    async function core(conn: AsyncDuckDBConnection, abort: AbortSignal) {

      
      const sql = filterMode == "sql" ? filter.sql : makeSqlQuery(filter, queryFromTable(filter))
      console.log(sql)
      abort.throwIfAborted()
      await ensureViews(conn, abort, sql)
      if (false) {
        conn.query(sql).then(t => t.schema.metadata)
      }

      abort.throwIfAborted()
      resultsPromise = conn.query(sql)
      resultsAgg = {}
      const resultTable = await resultsPromise
      resultsCount = resultTable.numRows
      results = Array.from(convertQueryResults(resultTable, parsePairingType(selectedPairing), 100));
      console.log({resultsPromise, results})
      const cols = resultTable.schema.names
      if (cols.includes("type") && cols.includes("res1") && cols.includes("res2")) {
        abort.throwIfAborted()
        resultsAgg.types = Object.fromEntries(
          Array.from(await conn.query(aggregateTypesQuery(sql))).map(x => [x.type, x.count])
        )
        console.log("types", resultsAgg.types)
        resultsAgg = resultsAgg
      }
      if (cols.includes("pdbid")) {
        abort.throwIfAborted()
        resultsAgg.pdbStructures = Object.fromEntries(
          Array.from(await conn.query(aggregatePdbCountQuery(sql))).map(x => [x.pdbid, x.count])
        )
        resultsAgg = resultsAgg
      }
      if (cols.some(x => String(x).startsWith("hb_"))) {
        abort.throwIfAborted()
        const x = await conn.query(aggregateBondParameters(sql, cols.map(c => String(c))))
        resultsAgg.bondStats = []
        for (let n of x.schema.names) {
          n = String(n)
          const m = n.match(/^hb_(\d+)_(.*)_([a-z0-9]+)$/)!
          if (m == null) continue
          const [ _, bond, param, stat ] = m
          resultsAgg.bondStats.push({ bond: Number(bond), param: param as any, stat: stat as any, value: Number([... x][0][n]) })
        }
        console.log("stats", resultsAgg.bondStats)
        resultsAgg = resultsAgg
      }
    }

    try {
      const conn = await connPromise
      await tableLoadPromise
      // await conn.cancelSent()
      return await updateResultsLock.withCancellableLock(abortSignal => core(conn, abortSignal))

    } catch (e) {
      resultsPromise = Promise.reject(e)
      console.error("Loading data failed:", e)
      throw e
    }
  }

  function getPairTypeList(selectedFamily) {
    const list = db.pairTypes.filter(p => p[0].toLowerCase() == selectedFamily.toLowerCase() || selectedFamily == null)

    list.sort((a, b) => a[0].toLowerCase().localeCompare(b[0].toLowerCase()) || a[1].localeCompare(b[1]))
    return list
  }

  // const imgDir = "http://[2a01:4f8:c2c:bb6c:4::8]:12345/"
  const imgDir = "https://pairs.exyi.cz/img"
  // const imgDir = base+"/img"
</script>

<Modal>

<div class="selector buttons has-addons is-centered are-small" style="margin-bottom: 0px">
  {#each db.pairFamilies as family}
    <button
      class="button"
      class:is-info={selectedFamily == family}
      class:is-selected={selectedFamily == family}
      on:click={() => {
        if (selectedFamily == family)
          selectedFamily = null
        else {
          selectedFamily = family
          selectedPairing = selectedPairing.replace(/^[^-]*-/, `${family}-`)
        }
      }}
    >{family}</button>
  {/each}
</div>

<div class="selector buttons has-addons is-centered are-small">
  {#each getPairTypeList(selectedFamily) as [family, bases]}
    <button
      class="button"
      class:is-success={selectedPairing.toLowerCase() == `${family}-${bases}`.toLowerCase()}
      class:is-selected={selectedPairing.toLowerCase() == `${family}-${bases}`.toLowerCase()}

      on:click={() => {selectedPairing = `${family}-${bases}`}}>{selectedFamily == null ? family + "-" : ""}{bases}</button>
  {/each}
</div>
<div class="filters">
  <FilterEditor bind:filter={filter}
    selectingFromTable={filter && queryFromTable(filter)}
    metadata={getMetadata(selectedPairing)}
    bind:mode={filterMode} />
</div>
{#await resultsPromise}
<div style="display:flex; flex-direction: row;">
  <div style="flex-grow: 1;"></div>
  <Spinner></Spinner>
  <div style="flex-grow: 1;"></div>
</div>
{:then result}
  {#if filterMode == "sql"}
    <div>
      {#each result.schema.fields as field}
        <span class="tag is-light" class:is-success={recommendedColumns.includes(field.name) || requiredColumns.includes(field.name)}><b>{field.name}</b>: {field.type}</span>
      {/each}
    </div>
    <div>
      {#each requiredColumns as c}
        {#if !result.schema.fields.find(f => f.name == c)}
          <span class="tag is-light is-danger" class:is-danger={true}><b>{c}: MISSING</b></span>
        {/if}
      {/each}
    </div>
    <div>{result.numRows} results</div>

    {#if requiredColumns.some(c => !result.schema.fields.some(f => f.name ==c))}
      <table class="table is-narrow is-striped is-fullwidth">
        <thead>
          <tr>
            {#each result.schema.fields as f}
              <th>{f.name}</th>
            {/each}
          </tr>
          <tr>
            {#each result.schema.fields as f}
              <th>{f.type}</th>
            {/each}
          </tr>
        </thead>
        <tbody>
          {#each [...result] as row}
            <tr>
              {#each result.schema.fields as f}
                <td>{row[f.name]}</td>
              {/each}
            </tr>
          {/each}
      </table> 
    {/if}
  {/if}
{:catch error}
  <pre style="color: darkred">{error}</pre>
{/await}
{#if Object.keys(resultsAgg.types ?? 0)?.length}
  <div class="stats-row">
    {#each Object.entries(resultsAgg.types) as [type, count], i}
      {#if i > 0}, {/if}
      <strong>{count}</strong> × {type}
    {/each}
    {#if Object.keys(resultsAgg.pdbStructures ?? 0)?.length}
      <span title={Object.entries(resultsAgg.pdbStructures).map(([pdb, count]) => `${count} × ${pdb}`).slice(0, 20).join(", ") + (Object.keys(resultsAgg.pdbStructures).length > 20 ? ", …" : "")}>
        from <strong>{Object.keys(resultsAgg.pdbStructures).length}</strong> PDB structures
      </span>
    {/if}
  </div>
{/if}
<PairImages pairs={results} rootImages={imgDir} imgAttachement=".png" videoAttachement=".mp4" />
{#if resultsCount != results?.length && resultsCount > 0}
  <Spinner></Spinner>
  <p style="text-align: center">Loading more... than 100 items is not implemented at the moment</p>
{/if}

</Modal>
<style>
  .stats-row {
    margin-left: 1rem;
    margin-right: 1rem;
    border-top: 1px solid #ccc;
    border-bottom: 1px solid #ccc;
  }
/*
.selector .selected {
  border: 2px solid black;
}
.selector {

} */
</style>
