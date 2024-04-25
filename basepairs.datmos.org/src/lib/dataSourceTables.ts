import type { AsyncDuckDBConnection } from "@duckdb/duckdb-wasm"
import { type PairingFamily, tryParsePairingType, normalizePairType } from "$lib/pairing"
import metadata from "$lib/metadata"
import * as filterLoader from "$lib/predefinedFilterLoader"
import { makeSqlQuery } from "./dbModels"

export async function ensureViews(conn: AsyncDuckDBConnection, abort: AbortSignal, queryTables: Iterable<string>, selectedPairing: string) {
    let pairid
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
    if (tableSet.has('selectedpair_allcontacts_f') || tableSet.has('selectedpair_allcontacts_boundaries_f')) {
      await addView(selectedNorm, 'selectedpair_allcontacts_f', `${selectedNorm}-filtered-allcontacts`)
      tableSet.delete('selectedpair_allcontacts_f')
    }
    if (tableSet.has('selectedpair_allcontacts_boundaries_f')) {
      const f = filterLoader.addHBondLengthLimits(selectedNorm, filterLoader.toNtFilter(await filterLoader.defaultFilterLimits.value, selectedNorm, null))
      console.log("boundaries filter: ", f)
      await conn.query(`CREATE OR REPLACE VIEW 'selectedpair_allcontacts_boundaries_f' AS ${makeSqlQuery(f, 'selectedpair_allcontacts_f')}`)
      tableSet.delete('selectedpair_allcontacts_boundaries_f')
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
