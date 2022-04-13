Genetic code prediction from karyorelict and heterotrich ciliates
=================================================================

Among the ciliates, there is an unusually diverse number of genetic codes used
by different species, compared with other groups of eukaryotes. The
karyorelicts and heterotrichs have some of the most unusual types of codes,
where stop codons are potentially ambiguous and can also code for amino acids.

Aims:
 * Assemble single cell transcriptome RNAseq data from uncultivated ciliates
 * Assemble whole-genome-amplification DNAseq data from uncultivated ciliates
 * Predict genetic codes from above assembles
 * Annotate genes from above assemblies, identify clear examples of ambiguous
   stop genetic codes

See our [preprint
(doi:10.1101/2022.04.12.488043)](https://doi.org/10.1101/2022.04.12.488043) for
more information.


Data
----

 * Single-cell transcriptome RNAseq from uncultivated marine ciliates, collected from Roscoff, France
 * Published datasets downloaded from SRA


Data deposition
---------------

 * Primary read data from this study: [ENA PRJEB50648](https://www.ebi.ac.uk/ena/browser/view/PRJEB50648)
 * Snakemake pipelines:
   * [Transcriptome assembly (this repo)](https://github.com/Swart-lab/karyocode-workflow)
   * [BUSCO marker workflow](https://github.com/Swart-lab/karyocode-analysis-busco)
   * [PORC genetic codes workflow](https://github.com/Swart-lab/karyocode-analysis-porc)
   * [PORC software](https://github.com/Swart-lab/PORC)
