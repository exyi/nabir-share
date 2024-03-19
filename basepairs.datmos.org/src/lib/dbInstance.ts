
import { AsyncDuckDB, AsyncDuckDBConnection, DuckDBDataProtocol } from '@duckdb/duckdb-wasm';
import metadata from '$lib/metadata.ts';
import { initDB } from '$lib/duckdb'
import { normalizePairFamily, compareFamilies } from '$lib/pairing'
// let db: any = null;
let conn: AsyncDuckDBConnection | null = null

export const parquetFiles = {
}
export const pairTypes: [ string, string ][] = metadata.map(m => m.pair_type).filter(p => !p[0].startsWith('n'))
export const pairFamilies: string[] = [...new Set(pairTypes.map(t => normalizePairFamily(t[0])))]

pairFamilies.sort(compareFamilies)
const cacheBuster = '?v=8'

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
const defaultHost = "https://pairs.exyi.cz"
export const fileBase = (new URL('tables/', host == 'localhost' ? defaultHost : document.baseURI)).href
export const imgDir = host == 'localhost' ? `${defaultHost}/img` : (new URL('pregen-img', document.baseURI)).href

export function getConnectionSync(): AsyncDuckDBConnection {
  if (!conn)
    throw new Error("DuckDB connection not initialized")
  return conn
}

export async function connect(): Promise<AsyncDuckDBConnection> {
  if (conn) {
    return conn
  }
  console.log("LOADING DB")
  const db = await initDB();
  for (const [name, url] of Object.entries(parquetFiles)) {
    await db.registerFileURL(name, `${fileBase}${url}`, DuckDBDataProtocol.HTTP, false);
  }
  conn = await db.connect();

  // console.log("PREPARING VIEWS")
  // const existingTables = await db.getTableNames(conn, )
  // for (const [name, url] of Object.entries(parquetFiles)) {
  //   await conn.query(`CREATE OR REPLACE VIEW '${name}' AS SELECT * FROM parquet_scan('${name}')`);
  // }
  window["duckdbconn"] = conn
  return conn;
}
