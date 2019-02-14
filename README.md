# Stoichiometry drives eco-evolutionary feedbacks
This repository contains open-source code, data, and text files for the eco-evo-stoich project to evaluate the effects of nutrient stoichiometry on the eco-evolutionary feedbacks between a cyanobacterial host and its phage.

## Project Questions

## Repository contents
### data
**Population Dynamics (pop-dynamics)**

`ts-counts_raw`: epi-fluorescent microscope field densities for each of ten images collected for each chemostat through 100 samplings (s1-s100). Columns 1-6 correspond to information about the image, organism, chemostat, nutrient treatment replicate, phage treatment type (phage amended [Infect] or non-phage amended [Control], and nutrient treatment type (N-limited [N] or P-Limited[P])

`ts-counts`: reduced microscopy data (ts-counts_raw) using the function `microsope.counts`. Each value represents the mean of ten images for each sampling date. Columns 1-6 correspond to information about the microbe (*Synechococcus* or phage), unique chemostat identifier, treatment replicate, phage treatment type (phage amended [Infect] or non-phage amended [Control], and nutrient treatment type (N-limited [N] or P-Limited[P]).

`cid-means`: *Synechococcus* and phage population counts as long form

`cstat-means`: *Synechococcus* and phage population counts as matrix with calculated standard error of the mean for each population. This file was used for the production of Figure 1 and Figure S2.

**Evolutionary dynamics (evo-dynamics)**
`inf-mat`: binary infection data (0 = infection/lack of *Synechococcus* growth, 1 = no infection/*Synechococcus* growth) collected for each *Synechococcus* strain (columns) and phage strain (rows) challenge. Strain information including the treatment (trt; C = Control, T = Treatment), nutrient limitation (lim; N = nitrogen-limited, P = phosphorus-limited), the unique chemostat identifier (cID; Ni = nitrogen-limited replicate, Pi = phosphorus-limited replicate), the isolation day in the timecourse (daynumber), the series time point (tm.pt), and isolate number (iso). Each strain was assigned a unique strain identifier using the series time point, cID, and isolate number. NAs correspond to challenges with conflicting results in the triplicate analysis.  

`inf-mat-processed`: `inf-mat` file transformed from a matrix to long form data. This file was used for the production of all figures related to the infection data.

`20150710-infmatrixNcoevo`, `20150710-infmatrixPcoevo`, `20150710-infmatrixNcontrols`, `20150710-infmatrixPcontrols`: Calculated coevolutionary dynamics between coevolving and noncoevolving strains within the chemostats. Values represent the proportion of successful infections between time group isolates within treatment and control chemostats. These data were used for the production of Figure 2.

***Supplemental**

`PercChange`: Percent change in growth rates from triplicate cultures using the nutrient limited media (med.base: NL = N-limited, PL = P-limited) from the chemostat experiment with additional N (+N) or P (+P).  This data corresponds to Figure S1.

## Contributors
Dr. Jay T. Lennon and Dr. Megan L. Larsen
Department of Biology, Indiana University
