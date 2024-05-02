## Basepair browser webapp

Running on: https://basepairs.datmos.org,
Sometimes running on test env: https://pairs.exyi.cz


### Setup

* Install node.js, npm and Python
* Run `npm install` (in this directory)
* To start a development server locally: `npm run dev` and open `http://localhost:1922/`
    - large assets are fetched from the deployed instance
* To build a deployable static website: `npm run build`
    - parquet file are expected under `tables/`, that is generated using `pair_distributions.py --reexport=partitioned`
    - optionally, generated images are expected under `pregen-img/`, which can be generated using `gen_contact_images.py`

### Configuration

Some configurable knobs are in the src/lib/config.ts file. It includes the asset paths, CSVs with parameter ranges and some defaults.
