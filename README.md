# Rescoring based on orthogonal information to improve protein identification in MS/MS
Repository contains scripts implementing methods described in the paper 
"Utility of RNA-seq and GPMDB protein observation frequency for improving the sensitivity of
protein identification by Tandem MS"

Link: http://www.ncbi.nlm.nih.gov/pubmed/25026199

## Usage

##### bayesPipeline.R
Takes as input a suppl file, a decoy groups file and a parsed protXML file (parsed by running Abacus on verbose mode; http://abacustpp.sourceforge.net/) and outputs a file with adjusted probabilities for the protein identifications.

##### decoyGroups.tsv
Example decoyGroups.tsv file. The above script divides the decoys into 2 groups, using 1 group for estimating decoyProb
and the other group for estimating False Discovery Rate at the end. So the decoy groups file is just a quick way
to keep this division consistent between runs. This can also be replaced by doing some sort of consistently seeded
random sampling of the decoys. 

##### sample.suppl.tsv
Example suppl file. Basically contains for each protein id the orthogonal information, GPMDB prot identification frequency
& RNAseq RPKM, used for performing the re-scoring.

##### bayesPipeline_rswtd.R
Same as the bayesPipeline.R script, except that the sampling of target values to assign to decoys is weighted by probability
of target identification (higher the probability, lower the chance of it getting sampled). Probably more accurate but might
introduce bias. Results from this alternative pipeline are discussed in supplementary information of the paper. 

##### bayesPipeline_hyrs.R
Same as the bayesPipeline.R script, except that this uses hyperscore as the basis of ranking identifications and re-scoring
instead of iniprob (peptide probability) as used in the main script. Results from this pipeline provide proof of how the 
improvement is likely dependent on the grey area in the identifications. 

### Scripts for creating suppl file

##### supplFilesScripts/calcrpkm.R
Script computes RPKM for all transcripts in a given accepted_hits.bam file output from Tophat, using Bioconductor GenomicFeatures package. Requires an already created transcript DB or can connect to biomart to download a transcript db [comment or uncomment appropriate lines to change behavior].

##### supplFilesScripts/enst2ensp.R
Script maps transcript ids to protein ids, right now only for ensembl.

##### supplFilesScripts/gtf2txdb.R
Script can create a txdb (to use with calcrpkm script) from a GTF in case the annotations you want to use for transcripts are not present online.

##### supplFilesScripts/makeSuppl.R
Combines the ensp.rpkm.tsv file (from calcrpkm and enst2ensp scripts), a protlen file and the GPMDB protein frequencies file to create a suppl file. Note: the bayesPipeline script by default performs both RNAseq based and GPMDB based probability adjustments. But if only interested in using one type of orthogonal data modify the script and the suppl files accordingly.

##### supplFilesScripts/gpmdbHumanProteomeGuide.Apr30th2011.tsv
Sample GPMDB protein frequencies file. Generated on April 2011. To create a more up to date version, use files provided at ftp://ftp.thegpm.org/projects/annotation/human_proteome_guide/
