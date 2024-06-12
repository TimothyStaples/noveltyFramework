Code used to generate results, figures and tables for GEB submission: A conceptual framework for measuring ecological novelty.

The data used to generate this publication were obtained from the Neotoma Paleoecological Database: details are located in dataProcessing.R.

REPOSITORY TREE:

**./dataProcessing.R**: commented code to access Neotoma data and harmonize taxonomy. Raw outputs have been included in repository so this script does not need to be re-run to reproduce results.

**./novelShowcase.R**: commented code that reproduces analyses, figures and tables found in the main text.

**./functions**: Small ease-of-use functions read into top-levels scripts.

**./rawdata**: Input data objects.
* 	**./rawdata/Neotoma vascular plant records.rds**: storage of raw Neotoma download. See Neotoma docs for metadata.
* 	**./rawdata/processedRecords.rds**: Neotoma data with harmonized taxonomy. See Neotoma docs for metadata.
* 	**./rawdata/WFOTaxaTable.csv**: Reference table for taxa (unique "variablenames" in Neotoma records), with abundance ("count") and World Flora harmonized taxonomy for family, genus and species.
* 	**./rawdata/WFOTaxaTableEDITED.csv**: Version with manual adjustment of some taxa where World Flora could not parse records.
* 	**./rawdata/shapeFiles/**: USA shape files obtained from cenus.gov. Used for study region maps.
	
**./plots**: Figure outputs from novelShowcase.R

**./outputs**: Tabular outputs from novelShowcase.R. 

**./outputs/novelResults.csv**: Raw data used for figures.
* 	*siteid:* Neotoma site ID
* 	*long:* Longitude (decimal degrees)
* 	*lat:* Latitude (decimal degrees)
* 	*noAnalog, timeArrow, pastComp, presentComp:* Novelty values as per main text descriptions.
* 	"Which" columns: siteid for which reference set entity was used to define novelty measurements. Only used to link 	nMDS points in supplementary figure.
* 	*pastSamples, presentSamples:* Number of pollen samples averaged into the 500 year bin from 1500-2000AD 	("pastSamples") and 1000-1500AD ("presentSamples")
* 	*fig3Number:* Reference number for labelled sites in main text figures.