# Incidence of Pelvic Inflammatory Disease Associated with _Mycoplasma genitalium_ Infection: Evidence Synthesis of Cohort Study Data

_Joanna Lewis, Paddy J Horner and Peter J. White_

_National Institute for Health Research Health Protection Research Unit in Modelling Methodology and Medical Research Council Centre for Global Infectious Disease Analysis, Imperial College London School of Public Health, London, UK_
_Population, Policy and Practice, UCL Great Ormond Street Institute of Child Health, 30 Guilford Street, London, UK_
_National Institute for Health Research Health Protection Research Unit in Evaluation of Interventions, University of Bristol, UK_
_Population Health Sciences, University of Bristol, UK_
_Modelling and Economics Unit, National Infection Service, Public Health England, London, UK_

This repository contains R and Stan code to estimate the risk of pelvic inflammatory disease associated with chlamydia and _M. genitalium_ infection, in women enrolled in the Prevention of Pelvic Inflammation (POPI) trial (Oakeshott et al. _BMJ_ **340**, 2010; Oakeshott et al. _Clin Infect Dis_ **51**, 2010).

* `SIS-5.stan` contains the Stan model used for inference.
* `run-SIS-5.R` is an R script to run the Stan model and produce graphs and summaries of theoutput.
* `SIS-5.rds` contains sampled posterior distributions produced by running the Stan model.
