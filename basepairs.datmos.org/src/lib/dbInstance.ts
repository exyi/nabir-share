
import { AsyncDuckDB, AsyncDuckDBConnection, DuckDBDataProtocol } from '@duckdb/duckdb-wasm';
import metadata from '$lib/metadata.ts';
import { initDB } from '$lib/duckdb'
import { normalizePairFamily, compareFamilies } from '$lib/pairing'
// let db: any = null;
let conn: AsyncDuckDBConnection | null = null

export const parquetFiles = {
}
export const pairTypes = metadata.map(m => m.pair_type)
export const pairFamilies = [...new Set(pairTypes.map(t => normalizePairFamily(t[0])))]

pairFamilies.sort(compareFamilies)

for (const [family, bases] of pairTypes) {
  parquetFiles[`${family}-${bases}`] = `${family}-${bases}.parquet`
  parquetFiles[`${normalizePairFamily(family)}-${bases}`] = `${family}-${bases}.parquet`
  parquetFiles[`${family}-${bases}-filtered`] = `${family}-${bases}-filtered.parquet`
  parquetFiles[`${normalizePairFamily(family)}-${bases}-filtered`] = `${family}-${bases}-filtered.parquet`
  parquetFiles[`n${family}-${bases}`] = `n${family}-${bases}.parquet`
  parquetFiles[`${family}-${bases}_n`] = `n${family}-${bases}.parquet`
}
export const host = window.location.hostname.match(/(^|[.])localhost$/) ? 'localhost' : window.location.hostname
export const fileBase = `${host == 'localhost' ? 'https://pairs.exyi.cz/' : document.baseURI}tables/`

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
