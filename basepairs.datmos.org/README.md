## Basepair browser webapp

Running on: https://basepairs.datmos.org
Sometimes running on test env: https://pairs.exyi.cz


### Setup

* Install node.js and npm and Python
* Run `npm install`
* To run dev version locally: `npm run dev` and open `http://localhost:1922/`
    - large assets are fetched from the deployed instance
* To build a deployable static website: `npm run build`
    - parquet file are expected under `tables/`, that is generated using `pair_distributions.py --reexport=partitioned`
    - optionally, generated images are expected under `img/`, which can be generated using `gen_contact_images.py`
