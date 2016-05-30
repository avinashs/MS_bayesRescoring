#Feb 11th 2013, Avinash Shanmugam
#Script will combine a ensp.rpkm.tsv file with the gpmdb-protfreq file
#and protlen file to make a suppl file

args = commandArgs(trailingOnly = TRUE);

if(length(args) < 2)
{
	print("USAGE: Rscript makeSuppl.R <ensp.rpkm file> <outfile>");
	stop();
}

rpkmFile = args[1];
outfile = args[2];

gpmdbFile = "~/vcap.rnaseq/input.files/gpmdbHumanProteomeGuide/gpmdbHumanProteomeGuide.Apr30th2011.tsv";

protlenFile = "~/vcap.rnaseq/input.files/ensp66.ids.protlen.tsv";

#Read in rpkmfile

rpkm = read.table(rpkmFile,sep="\t",header=T,as.is=T);

names(rpkm) = c("protid","txid","rpkm");

gpm = read.table(gpmdbFile,sep="\t",header=T,as.is=);

names(gpm) = c("protid","gpmNobs","bestloge");

protlen = read.table(protlenFile,sep="\t",header=T,as.is=T);

names(protlen) = c("protid","protlen");

#Merge protlen and gpm dfs
suppl = merge(protlen,gpm,by="protid");

#Merge the rpkm df in too
suppl = merge(suppl,rpkm,by="protid");

write.table(suppl,file=outfile,sep="\t",row.names=F,quote=F);
