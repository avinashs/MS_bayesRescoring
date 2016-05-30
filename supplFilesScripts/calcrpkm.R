#Avinash Shanmugam,March 5th 2012
#This script will use Bioconductor packages (Rsamtools,GenomicFeatures)
#to read in a tophat aligned bam file and calculate the RPKM
#for all known Ensembl transcripts[from a genomicFeatures transcriptdb]

library(Rsamtools);
library(GenomicFeatures);

#File name of tophat o/p bam file
tph.file = "/data/projects/vcap.rnaseq/tophat.out/k562/accepted_hits.bam";

#Filename of transcript db sqlite file
txdb.file = "/data/projects/vcap.rnaseq/input.files/Ens66.hg19.txdb.sqlite";

#Filename of outfile to write rpkm dataframe
out.file = "/data/projects/vcap.rnaseq/suppl.files/k562/k562.enst.rpkm.tsv";

#Read in tophat alignment files
print("Reading in tophat alignments..");
tphAligns = readBamGappedAlignments(tph.file);

#Read in transcript db
print("Reading in transcript db..");
txdb = loadFeatures(txdb.file);

#In case there isn't already a downloaded transcript db sqlite
#we need to connect to ensembl and create it.
#Uncomment the following 3 lines & comment the line above to do this
#print("Creating Transcript db; This will take a while..!");
#txdb = makeTranscriptDbFromBiomart(biomart="ensembl",dataset="hsapiens_gene_ensembl",host="aug2014.archive.ensembl.org");
#print("Transcript db created!");
#saveFeatures(txdb,file=txdb.file);

#Extract out coordinates of the exons, grouped by transcript
txExons = exonsBy(txdb,"tx");

#Determine the no. of reads overlapping the exons
print("Finding overlap read counts..");
readCounts =countOverlaps(txExons,tphAligns,ignore.strand=TRUE);

#Delete tophat Alignments to free up memory
rm(tphAligns);

##Calculating RPKM
print("Calculating RPKM");

#First find the tr length in Kb
numBases = sum(width(txExons));
txLengthsInKb = numBases / 1000;

#Get total number of reads mapped in millions
millionsMapped = sum(readCounts) / 1000000;

readsPerMillion = readCounts / millionsMapped;

#Get rpkm
rpkm = readsPerMillion /txLengthsInKb;

print("Creating output file..");

#Extract the transcriptids corresponding to the rpkm
txs = transcripts(txdb, columns=c("tx_id","tx_name") );

#Note: The transcripts fn doesn't retrieve the transcripts in order
#of the transcript ids. So before binding transcript ids and rpkm vals
#we need to sort them back into correct order

txnames.df = elementMetadata(txs)[,c("tx_name","tx_id")];

sort.txnames = txnames.df[with(txnames.df,order(tx_id)),];

enstid = as.vector(unlist(sort.txnames$tx_name));

#Create a data frame to print out the result
rpkm.df = data.frame(enstid,rpkm,readCounts,txLengthsInKb);

#Write to outfile
write.table(rpkm.df,file=out.file,sep="\t",row.names=FALSE,quote=FALSE);

