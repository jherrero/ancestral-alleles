AncestralAlleles
================

This software is intended for the inference of Ancestral Alleles (AA) using a
whole-genome multiple alignment.

Usage
-----

perl run_aa.pl --emf EMF_DIR --input INPUT.txt --output OUTPUT.txt

This will go through all the variants in the INPUT.txt file and calculates the AA
for each of them. At the moment only a specific TSV file format is supported, but
new file format (notably VCF will be supported in the near future).

Dependencies
------------

- Ortheus: https://github.com/benedictpaten/ortheus, although
 https://github.com/jherrero/ortheus/tree/original-code is recommended
- Muscle: http://www.drive5.com/muscle/


Credits
-------
Much has chenged since the previous verisons, but this is largely inspired by the
work done by Kathryn Beal while at the European Bioinformatics Institute
(http://www.ebi.ac.uk). This was done under the supervision of Javier Herrero,
in Paul Flicek's group.

The original files can be found in the Ensembl Compara repository, under
modules/Bio/EnsEMBL/Compara/RunnableDB/AncestralAllelesForIndels 
