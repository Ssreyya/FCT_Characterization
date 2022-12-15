# FCT_Characterization

This code can be used to extract the FCT Locus from SDSE genomes using their Annotated and Assembly files. Alternatively, it can be modified slightly 
for extraction of any locus (provided that locus has uniform flanking regions in all isolates). This can be processed in a linux environment.

Inputs:
1. ID_SDSE.txt - A text file containing names of all samples used. (Example: SRRxxxx)
2. A directory "Assemblies" in which all Assembly files are present. (Example: SRRxxxx.fasta)
3. A directory "Annotated" in which all Annotated files are present. (Example: SRRxxxx.gbk)
4. Empty directories named FCTgenes and temp
5. Ref_Ann_genes.fasta containing all the reference genes.
6. Blast databases must be present (check the path for it) to perform blastx  for the genes not having hits with genes in Ref_Ann_genes.fasta


Outputs:
1. FCT locus genes multi fastafiles for all provided isolates are created in the FCTgenes directory (Example: SRRxxxFCTgenes.fasta)
2. coregenes.txt contains all the coregenes for each isolate in comma-delimited format.
3. genelistnew.txt contains all the genes for all isolates. 


