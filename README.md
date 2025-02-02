## README

# This repository contains data and code from a study of phenotypic plasticity in Sorex shrews.

# REFERENCES
Bryan S. McLean, Kristin E. Stierman, Leo R. Ivey, Amanda K. Weller, Olivia S. Chapman, Ava M. Miller, Abigail Mendoza Garcia, Jada S. Byrd, Stephen E. Greiman. "Seasonal body size plasticity and the generality of Dehnel's Phenomenon in Sorex shrews." _American Naturalist_.

# DATA SETS
_SOCI_metadata_crania_femora_ageclass.csv_ - Morphological trait data for S. cinereus from Yancey County, North Carolina, USA. Individual specimens are indexed by institution and catalog number, with globally unique identifiers (Darwin Core triplets). Field-derived measurements presented here include body mass (in grams); total body length and tail length (both in mm); and sex. Additional attributes for specimens are available on the Arctos collection database. Skull measurements (derived from CT scans) include braincase height (BCH), braincase width (BCW), condylobasal length (SKL), and molariform toothrow length (TR; all in mm). Femur measurements (derived from CT scans) include femur total length (TL), femur width (W), and femur head diameter (HD). Tooth wear class (a proxy for age) is also given and follows classifications of J. O. Whitaker, Sorex cinereus. Mamm. Species 743, 1–9 (2004). 

_data-extractions_and_metadata.csv_ - Body size and braincase height measurements plus observational metadata extracted from prior studies of Sorex. Data were extracted at seasonal scale according to the parameters described in text above. Also included are unique observational identifiers (assigned by the authors, linkable by author and year as compiled in the bibliography in this file), geocoordinates (from prior studies or conducted newly here), and taxonomic information.

# CODE (all written by B. S. McLean)
_dehnel-analysis.R_ - Code for statistical tests and plotting of seasonal trait data in Sorex cinereus, Yancey County, NC, USA.

_dehnel-metaanalysis.R_ - Code for bioclimate PCA, phylogenetic signal tests, and Bayesian modeling of body size and braincase height plasticity in Sorex from Europe, Asia, and North America.
