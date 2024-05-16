# plasma concentration conversion
The script was used in the data analysis of the described paper below. It uses the report.gg_matrix.tsv file from DIA-NN as an input file to convert the MaxLFQ values to g/ml concentrations. It's based on a linear regression model with 22 known reported average values of proteins in serum (A2M, B2M, C1R, C2, C6, C9, CFP, CP, F10, F12, F2, F7, F8, F9, HP, KLKB1, MB, MBL2, SERPINA1, TFRC, TTR, and VWF) extracted from *Schenk, S. et al. A high confidence, manually validated human blood plasma protein reference set. BMC Med Genomics 1, 41 (2008). https://doi.org/10.1186/1755-8794-1-41*

### Longitudinal Fluctuations in Protein Concentrations and Higher-Order Structures in the Plasma Proteome of Kidney Failure Patients Subjected to a Kidney Transplant

*Sofia Kalaidopoulou Nteak12, Franziska Völlmy12, Marie V Lukassen12, Henk van den Toorn12, Maurits A den Boer12, Albert Bondt12, Sjors P A van der Lans3, Pieter-Jan Haas3, Arjan D van Zuilen4, Suzan H M Rooijakkers3, Albert J R Heck12**

1. Biomolecular Mass Spectrometry and Proteomics, Bijvoet Center for Biomolecular Research and Utrecht Institute for Pharmaceutical Sciences, University of Utrecht, Utrecht 3584 CH, The Netherlands.
2. Netherlands Proteomics Center, Utrecht 3584 CH, The Netherlands.
3. Department of Medical Microbiology, University Medical Center Utrecht, Utrecht 3584 CH, The Netherlands.
4. Department of Nephrology and Hypertension, University Medical Center Utrecht, Utrecht University, Utrecht 3584 CH, The Netherlands.

*Corresponding author: Albert J. R. Heck, E-mail: a.j.r.heck@uu.nl
