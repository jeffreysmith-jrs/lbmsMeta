# data/

The analysis reads its inputs from this folder. None of them are stored in the repository.

1. **Database.** Download the dataset from Zenodo (DOI to be added on publication) and copy the files in
   its `data/` folder here:

   | File | Contents |
   |---|---|
   | `observations.csv` | One row per extracted comparison of species richness between an intervention and its reference habitat, all interventions |
   | `studies.csv` | Every study assessed at full text, with its outcome and citation |
   | `taxa_groups.csv` | Surveyed taxon to broad taxonomic group |
   | `paper_language.csv` | Language of every paper assessed at full text |
   | `title_abstract_language.csv` | Web of Science language of every record screened at title/abstract |
   | `data_dictionary.csv` | Description of every column |

2. **RESOLVE Ecoregions 2017** (Dinerstein, E. et al. 2017. An ecoregion-based approach to protecting
   half the terrestrial realm. *BioScience* 67:534-545). Download the shapefile from
   <https://ecoregions.appspot.com/> and unzip it so that the path is
   `data/Ecoregions2017/Ecoregions2017.shp`. It is used to assign each site's biome.
