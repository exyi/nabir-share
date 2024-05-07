import type { ComparisonMode, NucleotideFilterModel } from "./dbModels";

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

    parameterBoundariesUrls: {
        'googleTable': 'https://docs.google.com/spreadsheets/d/e/2PACX-1vTvEpcubhqyJoPTmL3wtq0677tdIRnkTghJcbPtflUdfvyzt4xovKJxBHvH2Y1VyaFSU5S2BZIimmSD/pub?gid=245758142&single=true&output=csv',
        fr3dData: 'filters/fr3d-data-boundaries.csv'
    },
    defaultBoundaries: [
        'googleTable',
        'fr3dData'
    ]
}
