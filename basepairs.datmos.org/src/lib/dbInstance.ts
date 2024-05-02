
import { AsyncDuckDB, AsyncDuckDBConnection, DuckDBDataProtocol } from '@duckdb/duckdb-wasm';
import metadata from '$lib/metadata';
import { initDB } from '$lib/duckdb'
import { normalizePairFamily, compareFamilies } from '$lib/pairing'
import config from '$lib/config';
// let db: any = null;
let conn: AsyncDuckDBConnection | null = null
let db: AsyncDuckDB | null = null

export const parquetFiles: { [s: string]: string } = {
}
export const pairTypes: [ string, string ][] = metadata.map(m => m.pair_type as [string, string]).filter(p => !p[0].startsWith('n'))
export const pairFamilies: string[] = [...new Set(pairTypes.map(t => normalizePairFamily(t[0])))]
const registeredFiles = new Set<string>()

pairFamilies.sort(compareFamilies)
const cacheBuster = '?v=11'

for (const pairMeta of metadata) {
  const [family, bases] = pairMeta.pair_type
  const nf = normalizePairFamily(family)
  if (pairMeta.count != 0) {
    // parquetFiles[`${family}-${bases}`] = `${nf}-${bases}.parquet${cacheBuster}`
    parquetFiles[`${nf}-${bases}`] = `${nf}-${bases}.parquet${cacheBuster}`
    // parquetFiles[`${family}-${bases}-filtered`] = `${nf}-${bases}-filtered.parquet${cacheBuster}`
    parquetFiles[`${nf}-${bases}-filtered`] = `${nf}-${bases}-filtered.parquet${cacheBuster}`
    parquetFiles[`${nf}-${bases}-filtered-allcontacts`] = `${nf}-${bases}-filtered-allcontacts.parquet${cacheBuster}`
    parquetFiles[`n${nf}-${bases}`] = `n${nf}-${bases}.parquet${cacheBuster}`
    parquetFiles[`${nf}-${bases}_n`] = `n${nf}-${bases}.parquet${cacheBuster}`
  }
}
export const host = window.location.hostname.match(/(^|[.])localhost$/) ? 'localhost' : window.location.hostname
export const assetBaseUri = host == 'localhost' ? config.debugHost : document.baseURI
export const fileBase = (new URL(config.tablesPath, assetBaseUri)).href
export const imgDir = (new URL(config.imgPath, assetBaseUri)).href

export function getConnectionSync(): AsyncDuckDBConnection {
  if (!conn)
    throw new Error("DuckDB connection not initialized")
  return conn
}

export async function ensureAllFiles() {
  await ensureFiles(Object.keys(parquetFiles))
}

export async function ensureFiles(files: string[]) {
  if (!db) {
    throw new Error("DuckDB connection not initialized")
  }

  for (const file of files) {
    if (!registeredFiles.has(file) && file in parquetFiles) {
      registeredFiles.add(file)
      const url = new URL(parquetFiles[file], fileBase.replace(/\/$/, '') + "/").href
      await db.registerFileURL(file, url, DuckDBDataProtocol.HTTP, false)
    }
  }
}

export async function connect(): Promise<AsyncDuckDBConnection> {
  if (conn) {
    return conn
  }
  console.log("LOADING DB")
  db = await initDB();
  await ensureAllFiles()
  conn = await db.connect();

  // console.log("PREPARING VIEWS")
  // const existingTables = await db.getTableNames(conn, )
  // for (const [name, url] of Object.entries(parquetFiles)) {
  //   await conn.query(`CREATE OR REPLACE VIEW '${name}' AS SELECT * FROM parquet_scan('${name}')`);
  // }
  window["duckdbconn"] = conn
  return conn;
}
