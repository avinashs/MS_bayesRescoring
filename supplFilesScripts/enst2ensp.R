#April 25th 2012, Avinash Shanmugam
#Script will take a tsv file containing enstids and corresponding values
#and using the biomaRt library, will translate enstids to enspids and
#write to an output file.
#Note: Assumes ensembl 66 human database as the biomart to use


enst2ensp = function(in.df)
{
	#Extract enstids alone as a vector
	enst = in.df$enstid;

	#Load biomaRt library
	library(biomaRt);

	#Connect to ensembl
	ensembl = useMart(host="feb2012.archive.ensembl.org",biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl");

	#Query the mappings from biomart
	print("Querying mappings from biomart");
	ensp.map = getBM(attributes=c("ensembl_peptide_id","ensembl_transcript_id"),filters="ensembl_transcript_id", values=enst, mart=ensembl);

	#Change names of the mapped dataframe
	names(ensp.map) = c("enspid","enstid");

	#Remove rows that have a blank enspid 
	#These correspond to non-translated transcripts(rna-genes?)
	ensp.map = ensp.map[ensp.map$enspid != "",];

	#Merge the map df with original dataframe
	print("Merging mappings with the input dataframe");
	out.df = merge(ensp.map,in.df,by="enstid");

	#During merge, the function puts the column we are merging by as
	#the first column. But we want the enspid column to appear first
	#So just switching the 1st and 2nd cols
	out.df = out.df[,c(2,1,3:length(names(out.df)))];

	return(out.df);
}

############################## Main ###################################

args = commandArgs(trailingOnly =T);

if(length(args) < 2)
{
	print("USAGE: Rscript enst2ensp.R infile.tsv outfile.tsv");
	stop();
}

#Parse and store in and outfile names in variables
infile = args[1];
outfile = args[2];

#Read in the infile
print(paste("Loading input file:",infile));
in.df = read.table(infile,sep="\t",header=TRUE,as.is=TRUE);

out.df = enst2ensp(in.df[,c(1,2)]);

#write the out df to outfile
print(paste("Writing to outfile..",outfile));
write.table(out.df,file=outfile,sep="\t",row.names=FALSE,quote=FALSE);
