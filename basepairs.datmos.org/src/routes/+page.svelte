<script lang="ts">
  import PairImages from '$lib/components/pairimages.svelte'
  import Spinner from '$lib/components/Spinner.svelte'
  import metadata from '$lib/metadata'
	import FilterEditor from '$lib/components/filterEditor.svelte';
	import { aggregatePdbCountQuery, aggregateTypesQuery, defaultFilter, filterToSqlCondition, makeSqlQuery, parseUrl, type NucleotideFilterModel, filterToUrl, type StatisticsSettingsModel, statPresets, type StatPanelSettingsModel, statsToUrl, getDataSourceTable, makeDifferentialSqlQuerySingleTable, makeDifferentialSqlQuery, type ComparisonMode } from '$lib/dbModels.js';
	import { parsePairingType, type NucleotideId, type PairId, type PairingInfo, type HydrogenBondInfo, tryParsePairingType, type PairingFamily, normalizePairType } from '$lib/pairing.js';
	import { Modal } from 'svelte-simple-modal';
	import type { AsyncDuckDBConnection } from '@duckdb/duckdb-wasm';
	import { AsyncLock } from '$lib/lock.js';
  import * as db from '$lib/dbInstance';
  import * as arrow from 'apache-arrow'
	import { AsyncDebouncer } from '$lib/debouncer.js';
  import HistogramPlot from '$lib/components/HistogramPlot.svelte';
	import StatsPanel from '$lib/components/StatsPanel.svelte';
  import _ from 'lodash'
	import OverviewTable from '$lib/components/OverviewTable.svelte';

  let selectedFamily: string | undefined
  let selectedPairing: string | undefined
  const totalRowLimit = 30_000

  let filterMode: "basic"|"ranges" | "sql" = "basic"
  let filter: NucleotideFilterModel = defaultFilter()
  let filterBaseline: NucleotideFilterModel | undefined = undefined
  let comparisonMode: ComparisonMode = "difference"
  let miniStats: StatPanelSettingsModel[] = [ statPresets.histL, statPresets.histDA, statPresets.histAA ]
  let testStats: StatisticsSettingsModel = {
    enabled: false,
    panels: _.cloneDeep([ statPresets.histL, statPresets.histDA, statPresets.histAA ])
  }

  let lastUrlUpdate = 0

  function updateUrlNow(opt: {alwaysReplace?:boolean}={}) {
    const params = filterToUrl(filter, filterBaseline, filterMode)
    statsToUrl(params, testStats)
    const url = selectedPairing == null ? '' : `${selectedPairing}/${params.toString()}`
    if (window.location.hash.replace(/^#/, '') != url) {
      if (!opt.alwaysReplace && lastUrlUpdate + 1000 < performance.now()) {
        history.pushState({}, "", '#' + url)
      } else {
        history.replaceState({}, "", '#' + url)
      }
      lastUrlUpdate = performance.now()
    }
  }
  function updateUrl() {
    // urlUpdateDebouncer.debounce(updateUrlNow)
    updateUrlNow()
  }
  function setModeFromUrl(url: string) {
    url = url.replace(/^#/, '')
    const x = parseUrl(url)
    filterMode = x.mode
    filter = x.filter
    filterBaseline = x.baselineFilter
    if (x.stats != null) {
      testStats = x.stats
    } else {
      testStats.enabled = false
    }
    if (x.pairType != null) {
      selectedFamily = db.pairFamilies.find(f => f.toLowerCase() == x.pairFamily.toLowerCase()) ?? x.pairFamily
      selectedPairing = db.pairTypes.map(([f, b]) => `${f}-${b}`).find(f => f.toLowerCase() == x.pairType.toLowerCase()) ?? x.pairType
    } else {
      selectedPairing = undefined
      selectedFamily = undefined
    }
    updateUrlNow({alwaysReplace: true})
  }

  setModeFromUrl(window.location.hash)
  window.addEventListener("hashchange", (ev: HashChangeEvent) => setModeFromUrl(window.location.hash))

  $: {
    selectedPairing, filterMode, filter, filterBaseline, testStats
    updateUrl()
  }

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
          updateResultsLock.abortRunning()
          await conn.cancelSent()
          await conn.query(`DROP VIEW IF EXISTS selectedpair`)
          await conn.query(`DROP VIEW IF EXISTS selectedpair_f`)
          await conn.query(`DROP VIEW IF EXISTS selectedpair_n`)
          await conn.query(`DROP VIEW IF EXISTS selectedpair_allcontacts_f`)
          await conn.query(`DROP VIEW IF EXISTS selectedpair_allcontacts`)
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

  async function ensureViews(conn: AsyncDuckDBConnection, abort: AbortSignal, queryTables: Iterable<string>) {
    async function addView(pairingType: string | [ PairingFamily, string ] | null | undefined, name: string, file: string) {
      pairingType = tryParsePairingType(pairingType)
      if (pairingType && !metadata.some(m => m.pair_type[0].toLowerCase() == pairingType[0].toLowerCase() && m.pair_type[1].toLowerCase() == pairingType[1].toLowerCase())) {
        throw new Error(`Pairing type ${pairingType.join('-')} is not defined.`)
      }
      console.log(`Loading ${file} into view ${name}`)
      abort.throwIfAborted()
      await conn.query(`CREATE OR REPLACE VIEW '${name}' AS SELECT (pdbid || '-' || model || '-' || chain1 || '_' || coalesce(alt1, '') || nr1 || coalesce(ins1, '') || '-' || chain2 || '_' || coalesce(alt2, '') || nr2 || coalesce(ins2, '')) as pairid, * FROM parquet_scan('${file}')`)
    }
    const tableSet = new Set(queryTables)
    const existingTables = new Set([...await conn.query("select view_name as name from duckdb_views() union select table_name as name from duckdb_tables()")].map(r => r.name))
    const selectedNorm = normalizePairType(selectedPairing)

    for (const e in existingTables) {
      tableSet.delete(e)
    }
    console.log("missing tables:", [...tableSet])
    if (tableSet.has('selectedpair')) {
      await addView(selectedNorm, 'selectedpair', `${selectedPairing}`)
      tableSet.delete('selectedpair')
    }
    if (tableSet.has('selectedpair_f')) {
      await addView(selectedNorm, 'selectedpair_f', `${selectedNorm}-filtered`)
      tableSet.delete('selectedpair_f')
    }
    if (tableSet.has('selectedpair_allcontacts_f')) {
      await addView(selectedNorm, 'selectedpair_allcontacts_f', `${selectedNorm}-filtered-allcontacts`)
      tableSet.delete('selectedpair_allcontacts_f')
    }
    if (tableSet.has('selectedpair_n')) {
      await addView(selectedNorm, 'selectedpair_n', `n${selectedNorm}`)
      tableSet.delete('selectedpair_n')
    }

    for (const t of tableSet) {
      let pair
      if ((pair = tryParsePairingType(t)) != null) {
        await addView(pair, t, t)
      } else if ((t.endsWith('_f') || t.endsWith('-f') || t.endsWith('_n') || t.endsWith('-n')) &&
        (pair = tryParsePairingType(t.slice(0, -2))) != null) {
        await addView(pair, t, (t.endsWith('n') ? 'n' : '') + normalizePairType(t.slice(0, -2)) + (t.endsWith('f') ? '-filtered' : ''))
      } else if (/[_-]allcontacts[_-]f/i.test(t) &&
        (pair = tryParsePairingType(t.slice(0, -'-allcontacts'.length)))) {
        await addView(pair, t, normalizePairType(t.slice(0, -'-allcontacts'.length)) + '-filtered-allcontacts')
      } else {
        if (!existingTables.has(t))
          console.warn(`Maybe missing table: ${t}?`)
      }
    }
  }

  type ResultsAggregates = {
    count?: number
    types?: { [key: string]: number },
    pdbStructures?: { [key: string]: number },
    bondStats?: { bond: number, stat: "nncount" | "mean" | "min" | "max" | "median" | "p10" | "p25" | "p75" | "p90" | "stddev", param: "length" | "acceptor_angle" | "donor_angle", value: number }[]
  }
  let resultsPromise: Promise<arrow.Table<any> | undefined> = new Promise(() => {})
  let results = []
  let resultsTable: arrow.Table<any> | null = null
  let resultsCount = 0
  let resultsAgg: ResultsAggregates = {}

  $: {
    filter, filterMode, filterBaseline, selectedPairing, comparisonMode
    updateResults()
  }

  const requiredColumns = [ "pdbid", "chain1", "nr1", "chain2", "nr2", ]
  const recommendedColumns = [ "model", "ins1", "alt1", "res1", "res2", "ins2", "alt2", "res2" ]

  function getMetadata(pairType) {
    if (pairType == null) return undefined
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

    let count = 0
    for (const r of rs) {
      if (limit != null && count >= limit)
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
          function numberMaybe(x) { return x ? Number(x) : null }
          hbonds.push({
            length: numberMaybe(length),
            donorAngle: numberMaybe(r[`hb_${i}_donor_angle`]),
            acceptorAngle: numberMaybe(r[`hb_${i}_acceptor_angle`]),
            OOPA1: numberMaybe(r[`hb_${i}_OOPA1`]),
            OOPA2: numberMaybe(r[`hb_${i}_OOPA2`]),
            label: meta?.labels[i],
          })
        }
      }
      const coplanarity = r.bogopropeller ?? r.coplanarity_angle
      const { comparison_in_baseline, comparison_in_current } = r
      const comparison = (comparison_in_baseline || false) != (comparison_in_current || false) ? Boolean(comparison_in_current || false) : undefined
      yield { id, hbonds, coplanarity, comparison, originalRow: r }
      count++
    }
  }

  const updateResultsLock = new AsyncLock()

  async function updateResults() {
    async function core(conn: AsyncDuckDBConnection, abort: AbortSignal) {
      let startTime = performance.now()
      const timing = { ensureViews: -1, query: -1, queryConvert: -1, aggPdbId: -1, aggTypes: -1, aggBondParams: -1 }
      
      const metadata = getMetadata(selectedPairing)
      const limit = testStats.enabled ? totalRowLimit : 3000
      let sql = filterMode == "sql" ? filter.sql : makeSqlQuery(filter, getDataSourceTable(filter), null)
      if (filterBaseline != null) {
        // const comparisonMode: ComparisonMode = "missing"
        if (filterMode != "sql" && filterBaseline.sql == null && getDataSourceTable(filter) == getDataSourceTable(filterBaseline)) {
          sql = makeDifferentialSqlQuerySingleTable(filter, filterBaseline, getDataSourceTable(filter), null, comparisonMode)
        } else {
          const baselineSql = filterBaseline.sql || makeSqlQuery(filterBaseline, getDataSourceTable(filterBaseline))
          sql = makeDifferentialSqlQuery(sql, baselineSql, null, filter.orderBy, comparisonMode)
        }
      }

      console.log(sql)
      abort.throwIfAborted()
      const queryTables = new Set(await conn.getTableNames(sql))
      const requiredSelectedPair = [...queryTables].some(t => t.toLowerCase().startsWith("selectedpair"))
      if (selectedPairing == null || requiredSelectedPair && (metadata == null || metadata.count == 0)) {
        console.warn("Pair type not defined:", selectedPairing)
        resultsPromise = Promise.resolve(undefined)
        results = []
        resultsTable = new arrow.Table()
        resultsCount = 0
        resultsAgg = {}
        return
      }
      abort.throwIfAborted();
      await ensureViews(conn, abort, queryTables)
      abort.throwIfAborted()
      timing.ensureViews = performance.now() - startTime
      startTime = performance.now()

      if (false) {
        conn.query(sql).then(t => t.schema.metadata)
      }

      abort.throwIfAborted()
      resultsTable = await conn.query(`SELECT * FROM (${sql}) LIMIT ${limit}`)
      abort.throwIfAborted()
      resultsAgg = {}
      timing.query = performance.now() - startTime
      startTime = performance.now()

      resultsCount = resultsTable.numRows
      results = Array.from(convertQueryResults(resultsTable, parsePairingType(selectedPairing), 100));
      timing.queryConvert = performance.now() - startTime
      startTime = performance.now()

      // console.log({resultsPromise, results})
      const cols = resultsTable.schema.names
      if (cols.includes("type") && cols.includes("res1") && cols.includes("res2")) {
        abort.throwIfAborted()
        resultsAgg.types = Object.fromEntries(
          Array.from(await conn.query(aggregateTypesQuery(sql))).map(x => [x.type, Number(x.count)])
        )
        resultsAgg.count = Object.values(resultsAgg.types).reduce((a, b) => a + b, 0)
        console.log("types", resultsAgg.types)
        resultsAgg = resultsAgg
        timing.aggTypes = performance.now() - startTime
        startTime = performance.now()
      }
      if (cols.includes("pdbid")) {
        abort.throwIfAborted()
        resultsAgg.pdbStructures = Object.fromEntries(
          Array.from(await conn.query(aggregatePdbCountQuery(sql))).map(x => [x.pdbid, x.count])
        )
        resultsAgg = resultsAgg
        timing.aggPdbId = performance.now() - startTime
        startTime = performance.now()
      }
      if (cols.some(x => String(x).startsWith("hb_"))) {
        abort.throwIfAborted()
        // const x = await conn.query(aggregateBondParameters(sql, cols.map(c => String(c))))
        resultsAgg.bondStats = []
        // for (let n of x.schema.names) {
        //   n = String(n)
        //   const m = n.match(/^hb_(\d+)_(.*)_([a-z0-9]+)$/)!
        //   if (m == null) continue
        //   const [ _, bond, param, stat ] = m
        //   resultsAgg.bondStats.push({ bond: Number(bond), param: param as any, stat: stat as any, value: Number([... x][0][n]) })
        // }
        // console.log("stats", resultsAgg.bondStats)
        resultsAgg = resultsAgg
        timing.aggBondParams = performance.now() - startTime
        startTime = performance.now()
      }

      console.log("Loading finished: ", timing)

      return resultsTable
    }

    try {
      const conn = await connPromise
      await tableLoadPromise
      console.log(`Running new query`)
      updateResultsLock.abortRunning()
      await conn.cancelSent()
      return await updateResultsLock.withCancellableLock(abortSignal => resultsPromise = core(conn, abortSignal))

    } catch (e) {
      resultsPromise = Promise.reject(e)
      const isOk = e.message?.includes("The operation was aborted") ?? false
      if (!isOk) {
        console.error("Loading data failed:", e)
        throw e
      }
    }
  }

  function getPairTypeList(selectedFamily: string | null) {
    const isSymetrical = !!selectedFamily && selectedFamily[1] == selectedFamily[2]
    const allPairTypes = new Set(
      selectedFamily == null ? [] :
      isSymetrical ? [ 'A-A', 'A-C', 'A-G', 'A-U', 'C-C', 'C-U', 'G-C', 'G-G', 'G-U', 'U-U' ]
                   : [ 'A-A', 'A-C', 'A-G', 'A-U', 'C-A', 'C-C', 'C-G', 'C-U', 'G-A', 'G-C', 'G-G', 'G-U', 'U-A', 'U-C', 'U-G', 'U-U' ]);
    const pairs =
      db.pairTypes.filter(p => selectedFamily == null || p[0].toLowerCase() == selectedFamily.toLowerCase())
                  .map(p => ({ family: p[0], bases: p[1], real: true, conventional: selectedFamily == null || allPairTypes.has(p[1].toUpperCase()) }))
    const existingPairs = new Set(pairs.map(p => p.bases.toUpperCase()))
    for (const p of allPairTypes) {
      if (!existingPairs.has(p.toUpperCase())) {
        pairs.push({ bases: p, family: selectedFamily, real: false, conventional: true })
      }
    }

    pairs.sort((a, b) => a.family.toLowerCase().localeCompare(b.family.toLowerCase()) || a.bases.localeCompare(b.bases))
    return pairs
  }

  function formatTitleForPair(family: string, bases: string, real: boolean, conventional: boolean) {
    if (!real) {
      return 'No examples of this basepair type were found'
    }
    const m = metadata.find(m => m.pair_type[0] == family && m.pair_type[1] == bases)
    if (m == null) return "Does not exist"
    return [`${m.pair_type.join(' ')}`,
            `Class ${m.bp_class}`,
            `${m.labels.length} H-bonds`,
            `${m.med_quality + m.high_quality} pairs`,
            conventional ? "" : `${bases} is NOT symmetrical to ${bases[2]}-${bases[0]}!`
      ].join(', ')
  }

  // const imgDir = "http://[2a01:4f8:c2c:bb6c:4::8]:12345/"
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
        if (selectedFamily == family) {
          selectedFamily = null
          selectedPairing = null
        } else {
          selectedFamily = family
          selectedPairing = selectedPairing?.replace(/^[^-]*-/, `${family}-`)
        }
      }}
    ><b>{family}</b></button>
  {/each}
</div>

{#if selectedFamily != null}
<div class="selector buttons has-addons is-centered are-small" style="margin-bottom: 0px">
  {#each getPairTypeList(selectedFamily) as p}
  <!-- class:is-light={selectedPairing.toLowerCase() != `${p.family}-${p.bases}`.toLowerCase()} -->
    {@const m = metadata.find(m => m.pair_type[0] == p.family && m.pair_type[1] == p.bases)}
    <button
      class="button"
      class:is-light-warning={!p.conventional && selectedPairing?.toLowerCase() != `${p.family}-${p.bases}`.toLowerCase()}
      class:is-success={selectedPairing?.toLowerCase() == `${p.family}-${p.bases}`.toLowerCase()}
      class:is-selected={selectedPairing?.toLowerCase() == `${p.family}-${p.bases}`.toLowerCase()}
      class:is-static={!p.real}
      disabled={!p.real}
      title={formatTitleForPair(p.family, p.bases, p.real, p.conventional)}
      on:click={() => {
        if (selectedPairing?.toLowerCase() == `${p.family}-${p.bases}`.toLowerCase()) {
          selectedPairing = null
        } else {
          selectedPairing = `${p.family}-${p.bases}`
        }
      }}><b>{selectedFamily == null ? p.family + "-" : ""}{p.bases}</b>{#if m}&nbsp;({m.med_quality + m.high_quality}){/if}</button>
  {/each}
</div>
{#if selectedPairing && selectedPairing[1] == selectedPairing[2] && selectedPairing.slice(1) != 'SS'}
  <div class="buttons-help">The {selectedFamily} family is symmetrical, for example <b>C-A</b> is equivalent to <b>A-C</b></div>
{/if}
{/if}
<div style="margin-bottom: 1rem"></div>
{#if selectedPairing == null}
  <OverviewTable families={selectedFamily == null ? undefined : [selectedFamily]} />
{:else}
<div class="filters">
  <FilterEditor
    bind:filter={filter}
    bind:filterBaseline={filterBaseline}
    bind:comparisonMode={comparisonMode}
    selectingFromTable={filter && getDataSourceTable(filter)}
    metadata={getMetadata(selectedPairing)}
    bind:mode={filterMode} />
</div>
{#await resultsPromise}
<div style="display:flex; flex-direction: row;">
  <Spinner></Spinner>
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
    <div>{result?.numRows} results</div>

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

  <div style="position: absolute; width: 200px; text-align: center; margin-left: -100px; left:50%">
    {#if testStats.enabled}
    <a style="text-align: center;" href="javascript:;" on:click={e => testStats.enabled = false }>▲ collapse plots ▲</a>
    {:else}
    <a href="javascript:;" on:click={() => testStats.enabled = true}>▽ expand plots ▽</a>
    {/if}
  </div>
  <div>
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
  {#if !testStats.enabled}
    <!-- <div class="mini-stats">
      {#each fillStatsLegends({ panels: miniStats, enabled: true }, getMetadata(selectedPairing)).panels as stat}
        <HistogramPlot data={resultsTable} settings={stat} />
      {/each}
    </div> -->
  {/if}
</div>

{#if testStats.enabled}
  {#if resultsCount == totalRowLimit}
    <div class="notification is-warning">
      Only the first {resultsCount.toLocaleString("mfe")} rows are counted in the statistics
    </div>
  {/if}
  <StatsPanel data={resultsTable} bind:settings={testStats} metadata={getMetadata(selectedPairing)} />
{/if}

{/if}
{#await resultsPromise}
  <!-- <Spinner></Spinner> -->
{:then _} 
  {#if resultsCount == 0}
    {#if getMetadata(selectedPairing) == null}
      <h5 class="title is-5" style="text-align: center">Sorry, the base pair is not defined.</h5>
      <!-- TODO: show rotated bases to illustrate that it's probably impossible -->
    {:else if filterBaseline != null}
      <h5 class="title is-5" style="text-align: center">The compared datasets are equal.</h5>
      <p class="subtitle is-6" style="text-align: center">You are in comparison mode - only results missing in the current or the baseline dataset are shown</p>
      <p class="subtitle is-6" style="text-align: center">
        <button class="button is-warning" on:click={() => filterBaseline = undefined}>Reset baseline, exit comparison mode</button>
      </p>
    {:else}
      <h5 class="title is-5" style="text-align: center">Sorry, no results matched your query</h5>
      <p class="subtitle is-6" style="text-align: center">You can try to loosen the filters
        {#if filterMode != 'sql' && filter.filtered && (filter.datasource ?? "fred-f").startsWith("fred-")}
          (choose "entire PDB" datasource, for example)
        {/if}
      </p>
      {#if filterToSqlCondition(filter).filter(c => c!="jirka_approves").length > 0}
      <p style="text-align: center">
        <button class="button is-warning" on:click={() => filter = defaultFilter(filter?.datasource)}>Reset filter</button>
      </p>
      {/if}
    {/if}
  {/if}
{/await}
<PairImages pairs={results} rootImages={db.imgDir} imgAttachement={filter.rotX ? "-rotX.png" : ".png"} videoAttachement=".webm" videoOnHover={!!filter.rotX} />
{#if results && resultsCount != results.length && resultsCount > 0}
  <div style="margin: 100px; text-align: center">
    <!-- <p>Loading more... than 100 items is not implemented at the moment</p> -->
    <p style="text-align: center">This is the first {results.length} cases, total count is {resultsAgg.count ?? resultsCount}</p>
    <button type="button" on:click={() => {
      results = [...results, ...convertQueryResults(resultsTable.slice(results.length), parsePairingType(selectedPairing), 100) ]
    }}>
      Load {Math.min(100, resultsCount - results.length)} more examples
    </button>
  </div>
{/if}
{/if}

</Modal>
<style>
  .buttons-help {
    text-align: center;
    font-size: 0.75rem;
    margin-top: -0.5rem;
  }
  .stats-row {
    margin-left: 1rem;
    margin-right: 1rem;
    border-top: 1px solid #ccc;
    border-bottom: 1px solid #ccc;
  }
  .is-light-warning {
    color: rgba(0, 0, 0, 0.7);
    background-color: #fffaeb;
  }

/*
.selector .selected {
  border: 2px solid black;
}
.selector {

} */

.mini-stats {
  max-height: 100px;
  display: flex;
  flex-direction: row;
  justify-content: center;
}
</style>
