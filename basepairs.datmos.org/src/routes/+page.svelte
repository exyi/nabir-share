<script lang="ts">
  import PairImages from '$lib/components/pairimages.svelte'
  import Spinner from '$lib/components/Spinner.svelte'
  import metadata from '$lib/metadata'
	import FilterEditor from '$lib/components/filterEditor.svelte';
	import { aggregatePdbCountQuery, aggregateTypesQuery, defaultFilter, filterToSqlCondition, makeSqlQuery, parseUrl, type NucleotideFilterModel, filterToUrl, type StatisticsSettingsModel, statPresets, type StatPanelSettingsModel, statsToUrl, getDataSourceTable, makeDifferentialSqlQuerySingleTable, makeDifferentialSqlQuery, type ComparisonMode, type DetailModalViewModel, aggregateComparisoonTypeQuery, aggregateCountsAcrossGroups } from '$lib/dbModels';
	import { parsePairingType, type NucleotideId, type PairId, type PairingInfo, type HydrogenBondInfo, type PairingFamily, normalizePairType, getMetadata, convertQueryResults } from '$lib/pairing';
	import { Modal, type Context } from 'svelte-simple-modal';
	import type { AsyncDuckDBConnection } from '@duckdb/duckdb-wasm';
	import { AsyncLock } from '$lib/lock.js';
  import * as db from '$lib/dbInstance';
  import * as arrow from 'apache-arrow'
	import { AsyncDebouncer } from '$lib/debouncer';
  import HistogramPlot from '$lib/components/HistogramPlot.svelte';
	import StatsPanel from '$lib/components/StatsPanel.svelte';
  import _ from 'lodash'
	import OverviewTable from '$lib/components/OverviewTable.svelte';
	import DetailModal from '$lib/components/DetailModal.svelte';
	import { getContext } from 'svelte';
	import ContextHack from '$lib/components/ContextHack.svelte';
  import * as filterLoader from '$lib/predefinedFilterLoader'
  import config from '$lib/config'
	import { ensureViews } from '$lib/dataSourceTables';

  let selectedFamily: string | undefined // 'cWW' / 'tWW' / ... / 'tSS'
  let selectedPairing: string | undefined // 'cWW-A-A' / 'cWW-A-C' / ...
  const totalRowLimit = 30_000
  let normalTableLimit = 100

  let filterMode: "basic"|"ranges" | "sql" = "basic"
  let filter: NucleotideFilterModel = defaultFilter()
  let filterBaseline: NucleotideFilterModel | undefined = undefined
  let comparisonMode: ComparisonMode = config.defaultComparisonMode
  let miniStats: StatPanelSettingsModel[] = [ statPresets.histL, statPresets.histDA, statPresets.histAA ]
  let statistics: StatisticsSettingsModel = {
    enabled: false,
    panels: _.cloneDeep([ statPresets.histL, statPresets.histDA, statPresets.histAA ])
  }

  let modalContext: Context
  
  let detailModal : undefined | DetailModalViewModel = undefined

  let lastUrlUpdate = 0

  function getUrlParams(filter: NucleotideFilterModel, filterBaseline: NucleotideFilterModel | undefined, filterMode: string, testStats: StatisticsSettingsModel) {
    const params = filterToUrl(filter, filterBaseline, filterMode)
    statsToUrl(params, testStats)
    return params.toString()
  }

  function updateUrlNow(opt: {alwaysReplace?:boolean}={}) {
    const url = selectedPairing == null ? '' : `${selectedPairing}/${getUrlParams(filter, filterBaseline, filterMode, statistics)}`
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
      statistics = x.stats
    } else {
      statistics.enabled = false
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
    selectedPairing, filterMode, filter, filterBaseline, statistics
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
          // await conn.query(`DROP VIEW IF EXISTS selectedpair`)
          // await conn.query(`DROP VIEW IF EXISTS selectedpair_f`)
          // await conn.query(`DROP VIEW IF EXISTS selectedpair_n`)
          // await conn.query(`DROP VIEW IF EXISTS selectedpair_allcontacts_f`)
          // await conn.query(`DROP VIEW IF EXISTS selectedpair_allcontacts`)
          // await conn.query(`DROP VIEW IF EXISTS selectedpair_allcontacts_boundaries_f`)
          await conn.query(`DROP VIEW IF EXISTS selectedpair; DROP VIEW IF EXISTS selectedpair_f; DROP VIEW IF EXISTS selectedpair_n; DROP VIEW IF EXISTS selectedpair_allcontacts_f; DROP VIEW IF EXISTS selectedpair_allcontacts; DROP VIEW IF EXISTS selectedpair_allcontacts_boundaries_f;`)
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


  type ResultsAggregates = {
    count?: number
    types?: { [key: string]: number },
    pdbStructures?: { [key: string]: number },
    comparison?: { inBaseline: number, inCurrent: number, inBoth: number },
    bondStats?: { bond: number, stat: "nncount" | "mean" | "min" | "max" | "median" | "p10" | "p25" | "p75" | "p90" | "stddev", param: "length" | "acceptor_angle" | "donor_angle", value: number }[]
  }
  let resultsPromise: Promise<arrow.Table | undefined> = new Promise(() => {})
  let results = []
  let resultsTable: arrow.Table | null = null
  let resultsSchema: arrow.Schema<any> | null = resultsTable?.schema
  let resultsAgg: ResultsAggregates = {}

  $: {
    filter, filterMode, filterBaseline, selectedPairing, comparisonMode
    updateResults()
  }

  $: {
    statistics; statistics.enabled;
    statisticsChanged()
  }

  function statisticsChanged() {
    if (!statistics.enabled) return;
    if (resultsTable == null || resultsTable.numRows < Math.min(resultsAgg.count ?? Infinity, totalRowLimit)) {
      updateResults()
      return
    }
    if (statistics.panels.flatMap(p => p.variables).map(v => v.column).some(c => resultsSchema.names.includes(c) && !resultsTable.schema.names.includes(c))) {
      updateResults()
      return
    }
  }

  const requiredColumns = [ "pdbid", "chain1", "nr1", "chain2", "nr2", ]
  const recommendedColumns = [ "model", "ins1", "alt1", "res1", "res2", "ins2", "alt2", "res2" ]

  const updateResultsLock = new AsyncLock()

  async function updateResults(resultsMappingLimit = config.imgLimit) {
    async function core(conn: AsyncDuckDBConnection, abort: AbortSignal) {
      let startTime = performance.now()
      const timing = { ensureViews: -1, query: -1, querySchema: -1, queryConvert: -1, aggCounts: -1, aggBondParams: -1 }
      
      const metadata = getMetadata(selectedPairing)
      const limit = statistics.enabled ? totalRowLimit : normalTableLimit
      const basicColumns = ["pdbid", "model", "chain1", "res1", "nr1", "alt1", "ins1", "chain2", "res2", "nr2", "alt2", "ins2", "symmetry_operation1", "symmetry_operation2"]
      const comparisonColumns = filterBaseline == null ? [] : ["comparison_in_current", "comparison_in_baseline"]
      const statColumns = !statistics.enabled ? [] : [
        ...new Set(statistics.panels.flatMap(p => p.variables).map(v => v.column)),
      ]
      const projectionColumns = [...basicColumns, ...comparisonColumns, ...statColumns].filter(c => c).join(", ")


      let sql = filterMode == "sql" ? filter.sql : makeSqlQuery(filter, getDataSourceTable(filter), null)
      if (filterBaseline != null) {
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
        resultsSchema = resultsTable.schema
        resultsAgg = {}
        return
      }
      abort.throwIfAborted();
      await ensureViews(conn, abort, queryTables, selectedPairing)
      abort.throwIfAborted()
      timing.ensureViews = performance.now() - startTime; startTime = performance.now()

      abort.throwIfAborted()
      resultsTable = await conn.query(
        filterMode != "sql" ? `SELECT ${projectionColumns} FROM (${sql}) LIMIT ${limit}` :
        sql.trim().toLowerCase().startsWith('select') ? `SELECT * FROM (${sql}) LIMIT ${limit}` : sql
      )
      timing.query = performance.now() - startTime; startTime = performance.now()
      resultsSchema = resultsTable.schema
      abort.throwIfAborted()
      try {
        if (filterMode != "sql") {
          resultsSchema = (await conn.query(`SELECT * FROM (${sql}) LIMIT 0`)).schema
        }
      } catch { }
      timing.querySchema = performance.now() - startTime; startTime = performance.now()

      abort.throwIfAborted()
      resultsAgg = {}
      results = Array.from(convertQueryResults(resultsTable, parsePairingType(selectedPairing), resultsMappingLimit));
      timing.queryConvert = performance.now() - startTime; startTime = performance.now()

      // console.log({resultsPromise, results})
      const cols = resultsSchema.names
      const countGroupingSets: string[][] = [ [] ]
      if (cols.includes("family") && cols.includes("res1") && cols.includes("res2")) {
        countGroupingSets.push(["family", "ltrim(res1, 'D')", "ltrim(res2, 'D')"])
      }
      if (cols.includes("pdbid")) {
        countGroupingSets.push(["pdbid"])
      }
      if (cols.includes("comparison_in_baseline") && cols.includes("comparison_in_current")) {
        countGroupingSets.push(["comparison_in_baseline", "comparison_in_current"])
      }

      abort.throwIfAborted()
      const groupCountTable = await conn.query(aggregateCountsAcrossGroups(sql, countGroupingSets))
      for (const row of groupCountTable) {
        if (row.pdbid != null) {
          resultsAgg.pdbStructures ??= {}
          resultsAgg.pdbStructures[row.pdbid] = Number(row.count)
        }
        else if (row.family != null) {
          resultsAgg.types ??= {}
          resultsAgg.types[row.family + "-" + row["ltrim(res1, 'D')"] + "-" + row["ltrim(res2, 'D')"]] = Number(row.count)
        }
        else if (row.comparison_in_baseline != null) {
          resultsAgg.comparison ??= { inBaseline: 0, inCurrent: 0, inBoth: 0 }
          if (row.comparison_in_baseline && row.comparison_in_current) {
            resultsAgg.comparison.inBoth = Number(row.count)
          } else if (row.comparison_in_baseline) {
            resultsAgg.comparison.inBaseline = Number(row.count)
          } else if (row.comparison_in_current) {
            resultsAgg.comparison.inCurrent = Number(row.count)
          }
        }

        if (resultsAgg.types) {
          resultsAgg.count = Object.values(resultsAgg.types).reduce((a, b) => a + b, 0)
        } else if (resultsAgg.comparison) {
          resultsAgg.count = resultsAgg.comparison.inBaseline + resultsAgg.comparison.inCurrent + resultsAgg.comparison.inBoth
        } else if (resultsAgg.pdbStructures) {
          resultsAgg.count = Object.values(resultsAgg.pdbStructures).reduce((a, b) => a + b, 0)
        }
      }
      timing.aggCounts = performance.now() - startTime; startTime = performance.now()

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
    const isSymmetrical = !!selectedFamily && selectedFamily[1] == selectedFamily[2]
    const allPairTypes = new Set(
      selectedFamily == null ? [] :
      isSymmetrical ? [ 'A-A', 'A-C', 'A-G', 'A-U', 'C-C', 'C-U', 'G-C', 'G-G', 'G-U', 'U-U' ]
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

  function showDetailModal(d: DetailModalViewModel) {
    console.log("Opening modal", d)
    detailModal = d
    modalContext.open(DetailModal, {pair: d.pair, imageUrl: d.imgUrl, videoUrl: d.videoUrl, rotImageUrl: d.rotImgUrl, pairType: selectedPairing, filter, filterBaseline, requeryDB: filterMode != "sql" }, {
      classContent: "smodal-content",
      styleWindow: {
        width: "80vw",
      }
    })
  }

  // const imgDir = "http://[2a01:4f8:c2c:bb6c:4::8]:12345/"
  // const imgDir = base+"/img"
</script>

<Modal>

  <ContextHack name="simple-modal" bind:value={modalContext} />

<!-- <div style="float:left; padding: 0.5rem">
  <h1 class="title is-1"><a href="/">Home</a></h1>
</div> -->

<nav class="selector buttons has-addons is-centered" style="margin-bottom: 0px">
  <a class="button" href="/"
    on:click={() => {
      location.hash = ''
      window.location.reload()
    } }
    class:is-info={selectedFamily == null}
    class:is-selected={selectedFamily == null}
    style="color: {selectedFamily == null ? "white" : "red"}; font-weight: bold;"
    >Home</a>
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
</nav>

{#if selectedFamily != null}
<nav class="selector buttons has-addons is-centered" style="margin-bottom: 0px">
  {#each getPairTypeList(selectedFamily) as p}
  <!-- class:is-light={selectedPairing.toLowerCase() != `${p.family}-${p.bases}`.toLowerCase()} -->
    {@const m = metadata.find(m => m.pair_type[0] == p.family && m.pair_type[1] == p.bases)}
    <a
      class="button"
      href={`#${p.family}-${p.bases}/${getUrlParams(filter, filterBaseline, filterMode, statistics)}`}
      class:is-light-warning={!p.conventional && selectedPairing?.toLowerCase() != `${p.family}-${p.bases}`.toLowerCase()}
      class:is-success={selectedPairing?.toLowerCase() == `${p.family}-${p.bases}`.toLowerCase()}
      class:is-selected={selectedPairing?.toLowerCase() == `${p.family}-${p.bases}`.toLowerCase()}
      class:is-static={!p.real}
      class:is-disabled={!p.real}
      title={formatTitleForPair(p.family, p.bases, p.real, p.conventional)}
      on:click={ev => {
        if (selectedPairing?.toLowerCase() == `${p.family}-${p.bases}`.toLowerCase()) {
          selectedPairing = null
        } else {
          selectedPairing = `${p.family}-${p.bases}`
        }
        ev.preventDefault()
        return false
      }}><b>{selectedFamily == null ? p.family + "-" : ""}{p.bases}</b>{#if m && location.host != "xxbasepairs.datmos.org"}&nbsp;({m.med_quality + m.high_quality}){/if}</a>
  {/each}
</nav>
{#if selectedPairing && selectedPairing[1] == selectedPairing[2] && selectedPairing.slice(1, 3) != 'SS'}
  <div class="buttons-help">The {selectedFamily} family is symmetrical, <b>{selectedPairing.slice(4)}</b> is equivalent to <b>{selectedPairing.slice(4).split('-').toReversed().join('-')}</b></div>
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
          {#each [...result] as row, i}
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
    {#if statistics.enabled}
    <a style="text-align: center;" href="javascript:;" on:click={e => { statistics.enabled = false; return false } }>▲ collapse plots ▲</a>
    {:else}
    <a href="javascript:;" on:click={() => { statistics.enabled = true; return false }}>▽ expand plots ▽</a>
    {/if}
  </div>
  <div>
    {#if resultsAgg.comparison != null}
      {#if resultsAgg.comparison.inBaseline}
        <span class="tag is-light is-danger">only in baseline set: {resultsAgg.comparison.inBaseline}</span>
      {/if}
      {#if resultsAgg.comparison.inCurrent}
        <span class="tag is-light is-success">only in new set: {resultsAgg.comparison.inCurrent}</span>
      {/if}
      {#if resultsAgg.comparison.inBoth}
        <span class="tag is-light is-info">in both set: {resultsAgg.comparison.inBoth}</span>
      {/if}
      | 
    {/if}
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
  {#if !statistics.enabled}
    <!-- <div class="mini-stats">
      {#each fillStatsLegends({ panels: miniStats, enabled: true }, getMetadata(selectedPairing)).panels as stat}
        <HistogramPlot data={resultsTable} settings={stat} />
      {/each}
    </div> -->
  {/if}
</div>

{#if statistics.enabled}
  {#if resultsTable?.numRows == totalRowLimit}
    <div class="notification is-warning">
      Only the first {resultsTable.numRows.toLocaleString("mfe")} rows are counted in the statistics
    </div>
  {/if}
  <StatsPanel data={resultsTable} availableSchema={resultsSchema} bind:settings={statistics} metadata={getMetadata(selectedPairing)} />
{/if}

{/if}
{#await resultsPromise}
  <!-- <Spinner></Spinner> -->
{:then _} 
  {#if resultsTable?.numRows === 0}
    {#if getMetadata(selectedPairing) == null}
      <h5 class="title is-5" style="text-align: center">Sorry, the base pair is not defined.</h5>
      <!-- TODO: show rotated bases to illustrate that it's probably impossible -->
    {:else if filterBaseline != null}
      <h5 class="title is-5" style="text-align: center">
        {#if comparisonMode == "difference"}
          The compared datasets are equal.
        {:else if comparisonMode == "missing"}
          All entries from the current dataset are also present in baseline.
        {:else if comparisonMode == "new"}
          All entries in the baseline dataset are also present in the current dataset.
        {:else}
          Both datasets are empty.
        {/if}</h5>
      <p class="subtitle is-6" style="text-align: center">You are in comparison mode - only results missing in the current or the baseline dataset are shown</p>
      <p class="subtitle is-6" style="text-align: center">
        <button class="button is-warning" on:click={() => filterBaseline = undefined}>Reset baseline, exit comparison mode</button>
      </p>
    {:else}
      <h5 class="title is-5" style="text-align: center">Sorry, no results matched your query</h5>
      <p class="subtitle is-6" style="text-align: center">You can try to loosen the filters
        {#if filterMode != 'sql' && filter.filtered && (filter.datasource ?? config.defaultDataSource).startsWith("fr3d-")}
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
<PairImages pairs={results} rootImages={db.imgDir} rotImg={filter.rotX} imgAttachement=".png" videoAttachement=".webm" videoOnHover={!!filter.rotX}
  onClick={d => showDetailModal(d)} />
{#if results && resultsTable?.numRows > 0 && (resultsAgg.count ?? resultsTable.numRows) != results.length}
  <div style="margin: 100px; text-align: center">
    <p style="text-align: center">This is the first {results.length} cases, total count is {resultsAgg.count ?? resultsTable.numRows}</p>
    <button type="button" on:click={async () => {
      const oldCount = results.length
      const newCount = Math.min(config.imgLimit, (resultsAgg.count ?? resultsTable.numRows) - results.length)
      if (resultsTable.numRows < results.length + newCount) {
        normalTableLimit = results.length + newCount
        await updateResults(newCount + oldCount)
      }
      results = [...results, ...convertQueryResults(resultsTable.slice(results.length), parsePairingType(selectedPairing), oldCount + newCount - results.length) ]
    }}>
      Load {Math.min(config.imgLimit, (resultsAgg.count ?? resultsTable.numRows) - results.length)} more examples
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
