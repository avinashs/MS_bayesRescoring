#Avinash Shanmugam, March 22 2012
#Script will parse transcript information from a Genome annotation
#GTF file from Ensembl,put into dataframes of the necessary format,
#and then use makeTranscriptDb fn from Bioconductor:GenomicFeatures
#package to create a sqlite transcript db

gtf.file = "/data/projects/vcap.rnaseq/input.files/Homo_sapiens.GRCh37.62.gtf";
txdb.file = "/data/projects/vcap.rnaseq/input.files/Ens62.Hg19.txdb.sqlite";

##Read gtf file into a dataframe
print("Reading in gtf file..");
gtf.df = read.table(gtf.file,sep="\t",header=FALSE,as.is=TRUE);


#Store only selected cols into an intermediate dataframe
print("Creating a filtered intermediate dataframe..");
int.df = gtf.df[,c(1,3,4,5,7,9)];

names(int.df) = c("tx_chrom","feature","start","end","tx_strand","details");

#Filter out all rows that aren't in the chr list

chr.list =c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT");
int.df = int.df[int.df$tx_chrom %in% chr.list,];

#Extract ENST ids from the details string
print("Extracting ENST ids..");
tmp = as.character(sapply(int.df$details,function(x) unlist(strsplit(x,";"))[2]));

int.df$tx_name = as.character(sapply(tmp, function(x) unlist(strsplit(x," "))[3]));

##The first strsplit splits the details string on ";". We unlist the
##resulting list and take the 2nd element which is the transcript id.
##Now this is in the form " transcript_id ENST...". So we split this
##string on " " (spaces), unlist and take the 3rd element which is the 
##enstid string. (First element is just a trailing space).
##sapply is needed to apply this function separately for each row.
##as.character removes row names from the sapply result

#Extract Exon rank from the details string
print("Extracting exon rank..");
tmp = as.character(sapply(int.df$details,function(x) unlist(strsplit(x,";"))[3]));
int.df$exon_rank = as.character(sapply(tmp, function(x) unlist(strsplit(x," "))[3]));

#Extract ENSG ids from the details string
print("Extracting gene id..");
tmp = as.character(sapply(int.df$details,function(x) unlist(strsplit(x,";"))[1]));
int.df$gene_id = as.character(sapply(tmp, function(x) unlist(strsplit(x," "))[3]));


#Create transcript dataframe
print("Creating transcript dataframe..");
tx = data.frame(int.df$tx_name, int.df$tx_chrom, int.df$tx_strand);

tx = unique(tx);
names(tx)=c("tx_name","tx_chrom","tx_strand");

#Add tx_id column. tx_id is just a unique number for each transcript
tx$tx_id = c(1:length(tx$tx_name));

#Reorder columns so that tx_id is first column
tx = tx[,c(4,1,2,3)];

#Map back tx_id to int.df
int.df = merge(int.df, tx[,c(1,2)],by.x="tx_name",by.y="tx_name");


#Filling in tx_start and tx_end for the transcript df
#Create empty cols for the values first
tx$tx_start = 0;
tx$tx_end = 0;

for( i in c(1:length(tx$tx_name)))
{
	tmp.start = int.df$start[int.df$tx_id == i & int.df$feature =="exon"];
	tmp.end = int.df$end[int.df$tx_id ==i & int.df$feature =="exon"];

	#The start and end coordinates are already sorted in order of exon rank
	#in the gtf file. So to find the max or min, we only need to take the
	#first or last element in the vector.

	#For transcripts on + strand, the first exon's start and last exon's end
	#are the start and end of the transcript
	#But for transcripts on - strand, the exons are in rank order, while the
	#co-ordinates are still in + strand order. So beginning of transcript
	#is actually beginning of last exon and end is end of first exon

	if(tx$tx_strand[i] == "+")
	{
		tx$tx_start[i] = head(tmp.start,1);
		tx$tx_end[i] = tail(tmp.end,1);
	}
	else if(tx$tx_strand[i] == "-")
	{
		tx$tx_start[i] = tail(tmp.start,1);
		tx$tx_end[i] = head(tmp.end,1);
	}
}

##Creating exons dataframe
print("Creating exons dataframe..");
exons = int.df[int.df$feature == "exon",c(10,8,4,5)];
names(exons) = c("tx_id","exon_rank","exon_start","exon_end");

#It would also be useful to add CDS cols. But not sure how to handle exons for which there are no CDS (non-coding).
#So for now, I am just going to leave it as it is

##Creating genes dataframe
print("Creating genes dataframe..");
genes = merge(tx[,1],int.df[,c(9,10)],by.x="tx_id",by.y="tx_id");


print("Creating a txdb..");
library(GenomicFeatures);

txdb = makeTranscriptDb(tx, exons, genes);

saveFeatures(txdb,file=txdb.file);


