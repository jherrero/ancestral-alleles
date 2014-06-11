AncestralAlleles
================

This software is intended for the inference of Ancestral Alleles (AA) using a
whole-genome multiple alignment. The current focus in on 1-bp indels. In
addition, it can also take a FASTA file with pre-calculated ancestral sequence
to use for providing AA for SNPs.

Currently the repository only includes the software to generate the predictions
for all the possible indels at a given genomic position. The software requires
a pre-calculated whole-genome multiple alignment. It uses Ortheus for inferring
the ancestral sequences and Muscle to check for low-complex regions and filter
them out as the corresponding alignments would be unreliable.


Usage
-----

The inference module is designed to be used with the eHive code:

standaloneJob.pl RunAncestralAllelesOnEMF --emf Compara.6_primates_EPO.chr10_1.emf.gz \
 --vcf ALL.chr10.1bp_indels.vcf -anc homo_sapiens_ancestor_10.fa \
 --ortheus_exe ortheus_core --muscle muscle > Compara.6_primates_EPO.chr10_1.txt


Dependencies
------------

- eHive: https://github.com/Ensembl/ensembl-hive
- Ortheus: https://github.com/benedictpaten/ortheus
- Muscle: http://www.drive5.com/muscle/


Roadmap
-------

- [ ] VEP plugin using pre-calculated data
- [ ] Live VEP plugin, running the predictions on the fly 
- [ ] Support of MAF files (http://genome.ucsc.edu/FAQ/FAQformat.html#format5)


Credits
-------
This repository contains a re-write of the code initiated by Kathryn Beal while
at the European Bioinformatics Institute (http://www.ebi.ac.uk). This was done
under the supervision of Javier Herrero, in Paul Flicek's group.
