import type { ComparisonMode, NucleotideFilterModel } from "$lib/dbModels";

export default {
    defaultOrderBy: 'rmsdD',
    defaultDataSource: <NucleotideFilterModel["datasource"]>'allcontacts-f',
    defaultComparisonMode: <ComparisonMode>'union',
    /// default limit of displayed images in the gallery
    imgLimit: 100,
    /// hostname of server to use for images and parquet files when the webapp is running on localhost
    debugHost: 'https://pairs.exyi.cz',
    imgPath: '/pregen-img',
    tablesPath: '/tables',
}
