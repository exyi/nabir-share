import * as duckdb from '@duckdb/duckdb-wasm';
import duckdb_wasm from '@duckdb/duckdb-wasm/dist/duckdb-mvp.wasm?url';
import duckdb_worker from '@duckdb/duckdb-wasm/dist/duckdb-browser-mvp.worker.js?worker';

import type { AsyncDuckDB } from '@duckdb/duckdb-wasm';

let db: AsyncDuckDB | null = null;

const initDB = async () => {
	if (db) {
		return db; // Return existing database, if any
	}

	// Instantiate worker
	const logger = new duckdb.ConsoleLogger();
	const worker = new duckdb_worker();

	// and asynchronous database
	db = new duckdb.AsyncDuckDB(logger, worker);
	
	await db.instantiate(duckdb_wasm);
	await db.open({ query:  { castBigIntToDouble: true, castDecimalToDouble: true, }, filesystem: {allowFullHTTPReads: true } })
	return db;
};

export { initDB }; // so we can import this elsewhere
