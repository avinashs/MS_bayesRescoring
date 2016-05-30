#Avinash Shanmugam, Oct 17th 2012
#Functions and scripts to do a probWeighted sampling
#of forward values and assigning to decoys. This should give
#a more accurate distribution of values for the decoys


#Function will take an abacus out file with new.verbose output
#(outputs all the protIds, not just the representative proteins)
#and create extra cols to mark the representatives proteins for 
#every protId

markRepProts = function(df)
{
	#Non representative protein ids are specified in the form
	# RepProtId:::ProteinId
	#Make list of rows that have non representative protein ids
	nonRepRows = grep(":::",df$protid);


	#Make a column to store whether or not a particular protid
	#is a representative protein. Set default to 1.
	df$isrep = 1;

	#Set isrep=0 for all the nonRepRows
	df$isrep[nonRepRows] = 0;

	#Now split up the RepProtId:ProteinId strings into separate
	#repProtId and proteinid strings and store into a matrix
	splitIdstringsMatrix = matrix(unlist(strsplit(df$protid[nonRepRows],":::")),ncol=2,byrow=T);

	##strsplit returns a list with each element a vector with 2
	##elements, repProtId and protein Id. We unlist and convert it
	##into a matrix for ease of use.

	#Create a repProtid column. For the representative proteins
	#the repProtid will just be their protid itself. But for the
	#nonrep proteins, the protid of the repProtein as got from the
	#splitIdstrings. 

	df$repProtid = "";

	df$repProtid[-nonRepRows] = df$protid[-nonRepRows];
	df$repProtid[nonRepRows] = splitIdstringsMatrix[,1];

	#Also replace the repprotid:::proteinId string with just protid
	#for the nonRepRows
	df$protid[nonRepRows] = splitIdstringsMatrix[,2];

	return(df);
}

#Function will read in a decoyGroups file (file that marks the list 
#of all decoys into two groups). Using that file, it then marks decoys 
#from group -1, identified in df as -1 (isfwd value is set to -1).
markDprobDecoys = function(df,decoyGroupsFile)
{
	#Read in decoyGroupsFile
	dg = read.table(decoyGroupsFile,sep="\t",header=T,as.is=T);

	#Find the rows of decoys in df that have isfwd = -1 in dg
	dprobDecoys = df$protid %in% dg$protid[dg$isfwd ==-1];

	#Set isfwd for those decoys to -1 in df also
	df$isfwd[dprobDecoys] = -1;

	#Return df
	return(df);
}

#Function will calculate the false positive to decoy ratio (r-factor)
#Give an estimation of the number of forward proteins, 1 decoy is equivalent to

calcRfactor = function(df,decoyType,mxThreshold)
{
	#The maxIniProb threshold below which we look to estimate r
	mxThreshold = 0.2;

	#Get number of decoys and non-decoys (forwards + blinded decoys) below the mxThreshold
	ndecoys = nrow(df[df$maxiniprob <= mxThreshold & df$isrep ==1 & df$isfwd == decoyType,]);
	nforwards = nrow(df[df$maxiniprob <= mxThreshold & df$isrep ==1 & df$isfwd != decoyType,]);

	#Calculate r (false-positive to decoy ratio)
	r = nforwards / ndecoys;

	return(r);
}


#Function will calculate a decoy probability, in other words a localised
#False Discovery rate. This is the probability that any given protein
#is a false identification as estimated from the number of decoys in the
#maxIniprob bin in which that protein occurs.

calcDecoyProb = function(df, r)
{
	#No. of proteins to select in one bin when 
	#calculating local FDR
	binSize = 100;

	#Sort df in decreasing order of maxiniprob
	df = df[with(df,order(-maxhyperscore)),];
	
	#From df get subset of only the representative proteins
	repDf = df[df$isrep == 1,c("protid","isfwd","maxhyperscore")];
	
	#Create decoyProb and decoyProbRaw columns and set default to 1
	#(decoyProbRaw will store the original calculated vals and decoyProb
	#will store the value after applying loess smoothing)
	repDf$decoyProbRaw = 1;
	repDf$decoyProb = 1;

	#Get no. of rows in repDf
	dfEnd = nrow(repDf);

	#Set binStart and binEnd positions
	binStartPoints = seq(1,dfEnd,binSize);

	for(binStart in binStartPoints)
	{
		binEnd = binStart + binSize -1;

		#If binEnd is greater than dfEnd, i.e. there aren't enough
		#proteins remaining to make a full bin, just set binEnd to dfEnd
		if(binEnd > dfEnd)
		{
			binEnd = dfEnd;
		}

		#For the current bin, get the isfwd col alone from repDf
		isfwd = repDf$isfwd[c(binStart:binEnd)];

		#Get number of fwds and decoys
		nfwd = sum(isfwd !=-1);
		
		nrev = sum(isfwd ==-1);

		#Calculate the local fdr as 2*nrev/nfwd
		locFdr = (r*nrev)/nfwd;

		#Assign locFdr as decoyProb to all protein in current bin
		repDf$decoyProbRaw[c(binStart:binEnd)] = locFdr;
	}

	#Having calculated the raw decoy probabilities, we want to apply
	#Loess smoothing on the values to remove some of the localised 
	#fluctuations

	#But if there are a lot of decoys in the data, it will bias the
	#smoothing against high confidence ids. All ids with a decoy prob =0
	#might get a slightly higer value due to the smoothing. To prevent
	#this, we will only smooth the high / mid confidence region.
	#Low confidence ids will all be uniformly set to decoyProb = 1

	#The boundary of low confidence ids is assumed to be the bin at
	#which the decoyProbRaw reaches a value of 1 or higher
	if(max(repDf$decoyProbRaw) < 1)
	{
		smoothLimit = nrow(repDf);
	}
	else
	{
		#Find the point at which decoyProbRaw reaches 1 or above
		oneLimit = which(repDf$decoyProbRaw >=1)[1];

		#To get a continuous curve upon smoothing, include the rest of
		#the bin in which we reach 1 or above in the smoothLimit
		#Now, in cases where we pass dprob 1 at the last bin, we might not
		#have binSize number of proteins left over in the bin. 
		#So we set smooth limit to either oneLimit + binSize or the total no.
		#of proteins in repDf, whichever is lower
		smoothLimit = min((oneLimit + binSize -1),nrow(repDf));
	}

	#Perform loess smoothing
	loessModel = loess(decoyProbRaw~maxhyperscore, repDf[1:smoothLimit,]);

	#Use the loess model to predict smoothed dProb values
	dprobSmooth = predict(loessModel, repDf$maxhyperscore[1:smoothLimit]);

	#Create a col to store smoothed decoyProb values.Set default val 1
	repDf$decoyProb = 1;

	#Store the dprobSmooth values in the decoyProb col till smoothLimit
	repDf$decoyProb[1:smoothLimit] = dprobSmooth;

	#Now we have taken the localFdr to be the probability that
	#an identification from a particular bin is a decoy. But in
	#the low maxIniprob values, many bins contain more decoys than
	#fwds meaning localFdr > 1. So for decoyProb, we just set all
	#these values to 1. 

	#Doing this step after loess smoothing because Loess might
	#increase values that are 1 to greater than 1.
	repDf$decoyProb[repDf$decoyProb >1] = 1;
	repDf$decoyProb[repDf$decoyProb <0] = 0;


    #Now that the decoyProbs for all the repDf proteins have been
	#filled in, use this to fill in decoyProbs in the original Dfs
	#(Leave out isfwd, maxiniprob cols from repDf while merging 
	#,since those are already in the original df we don't want to duplicate them)

	df = merge(df, repDf[,-c(2,3)], by.x="repProtid", by.y="protid", all.x=TRUE,sort=FALSE);

	return(df);
}


addSupplFwdVals = function(df, sp)
{
	#Add the supplementary columns to df by merging df and sp dataframes
	df = merge(df,sp[,-2],by="protid",all.x=T);

	#Set any NA values in df to zero
	#Rev prots and contaminants.
	df[is.na(df)]=0;

	#Rev prot suppl values will get reassigned by later fns
	#But contaminant values will always remain 0

	return(df);
}

#Function will accept a suppl dataframe and an abacusOut dataframe with
#decoyProb calculated. It will merge and assign decoyProbs to proteins
#in suppl, with proteins absent in abacusOut being assigned decoyProb =1
#(i.e. Probability of being a false identification is 1)

addDprobToSuppl = function(sp,df)
{

	#Use merge to add decoyProb values to the suppl dataframe
	#Only selecting protid and decoyProb cols from df since we don't
	#need any of the other cols
	sp = merge(sp,df[,c("protid","decoyProb")],by="protid",all.x=T);

	#After the merge, proteins weren't identified in the mzXML file
	#will be assigned NA values. Set these values to 1 (i.e. definitely 
	#a false identification, decoyProb = 1)

	sp$decoyProb[is.na(sp$decoyProb)] = 1;

	return(sp);
}

#Function will take the df and sp dataframes, get the indexes of 
#sampling boundaries for length based sampling and return a boundsList.
#These bounds indexes can then be used with the assignSupplValsToRev fn
#to quickly assign values without having to iterate through and find
#the sampling boundaries each time

getSupplValsSamplingBounds = function(revProtlen,spProtlen)
{
	protlenRange = 50;
	minSampleSize = 100;

	#Get the longest protlen in df rev
	revMaxlen = tail(revProtlen,n=1);

	#Calculate what the last bin value must be based on maxlen
	revMaxbin = (trunc(revMaxlen/protlenRange) +1)*protlenRange;

	#Create bins to be used with hist
	revBins = seq(0,revMaxbin,by=protlenRange);

	#Get frequency counts in each bin
	#For sp counts, we will only consider protein that are less 
	#than or equal to the revMaxbin value (longest rev prot, rounded)
	revCounts = hist(revProtlen,breaks=revBins,plot=F)$counts;
	spCounts = hist(spProtlen[spProtlen <revMaxbin],breaks=revBins,plot=F)$counts;

	#We will be using the length of the revCounts vector in
	#several places following. So storing it in a variable so
	#that I don't have to call the length fn everytime.
	revCountsLength = length(revCounts);

	#Based on the frequency counts, we can get the indices
	#that will set the sampling boundaries for rev
	#First creating a blank aray with all zeros to store the indices
	spBounds = rep(0,revCountsLength);
	revBounds = rep(0,revCountsLength);

	boundsIndex = 1;
	countsIndex = 1;

	while(countsIndex <= revCountsLength)
	{
		grpEnd = countsIndex;
		grpSum = spCounts[countsIndex];

		while(grpSum < minSampleSize & grpEnd < revCountsLength)
		{
			grpEnd = grpEnd + 1;
			grpSum = sum(spCounts[countsIndex:grpEnd]);
		}
		
		if(grpSum >= minSampleSize)
		{
			spBounds[boundsIndex] = sum(spCounts[1:grpEnd]);
			revBounds[boundsIndex] = sum(revCounts[1:grpEnd]);
		}
		else
		{
			#If the grpSum is less than minSampleSize even after
			#going through the while loop, it means that there are no
			#proteins left to increase sampleSize. So we add all remaining
			#proteins into the previous group
			spBounds[boundsIndex -1] = sum(spCounts);
			revBounds[boundsIndex -1] = sum(revCounts);
		}

		boundsIndex = boundsIndex + 1;
		countsIndex = grpEnd + 1;
	}

	#Now that all the sampling boundaries have been added to the bounds vectors.
	#any remaining spaces on the vectors (which will have zero values) can 
	#be removed.
	spBounds = spBounds[-which(spBounds ==0)];
	revBounds = revBounds[-which(revBounds ==0)];

	#To make it easier to iterate through values in the next fn, we will add
	#zeros to the beginning of both these vectors
	spBounds = c(0,spBounds);
	revBounds = c(0,revBounds);

	boundsList = list(spBounds,revBounds);

	return(boundsList);

}

#Function uses the bounds list derived by getSupplValsSamplingBounds fn
#,randomly samples values from the sp dataframe according to the
#bounds and assigns them to decoys in the rev dataframe. The rev df with
#suppl vals added is returned.

assignSupplValsToRev = function(rev,sp,boundsList)
{
	#Col names of the suppl vals columns
	supplValsCols = c("gpmNobs","rpkm");

	#Sort the dataframes in increasing order of protlen
	rev=rev[order(rev$protlen),];
	sp = sp[order(sp$protlen),];

	#Extract out the vectors storing samplingBounds from the boundsList
	revBounds = boundsList[[2]];
	spBounds = boundsList[[1]];

	#Loop through the bondsVectors and perform the sampling
	for(i in c(2:length(revBounds)) )
	{
		revSampleIndices = c( (revBounds[[i-1]] + 1) : revBounds[[i]] );
		spSampleIndices = c( (spBounds[[i-1]] +1) : spBounds[[i]] );

		nRevSample = length(revSampleIndices);
		nSpSample = length(spSampleIndices);

		sampledRows = sample(spSampleIndices,nRevSample);

		rev[revSampleIndices,supplValsCols] = sp[sampledRows,supplValsCols];
	}

	return(rev);
}

#Function will take the bin thresholds supplied and each protein on
#which bin its logrpkm or loggpm value falls into

markBins = function(df,binningCol,binThresholds)
{
	#Creating the binCols colName
	binCol = paste(binningCol,"Bins",sep="");

	#Add a new col to df, fill it with zeros and name it binCol
	df$newCol = 0;
	names(df)[names(df) == "newCol"] = binCol;

	#Append -Inf and Inf to the binThresholds vector for ease of 
	#iterating through the thresholds
	binThresholds = append(-Inf, binThresholds);
	binThresholds = append(binThresholds, Inf);

	nBinThresholds = length(binThresholds);

	#Fill in bin values
	for( i in c(2:nBinThresholds))
	{	
		df[df[,binningCol] >= binThresholds[i -1] & df[,binningCol] < binThresholds[i], binCol] = i -1;
	
	}

	return(df);
}

#Function will take df with bin vals marked and computed adjP based
#on its decoy probability and conditional probabilities for its bin

calcAdjProb = function(df, binCol)
{
	#Get list of unique bin values
	binVals = unique(df[,binCol]);

	#Creating columns to store the probability given F
	#and probability given R values that will be calculated
	df$pgvnFcol=0;
	df$pgvnRcol=0;

	#Calculating number of forward and decoys in the data
	#We will only count forwards with decoyProb greater than 0.95
	#as forwards. While all decoys will be counted as decoys.
	nF = nrow(df[df$isfwd != -1 & df$decoyProb <= 0.5,]);
	nR = nrow(df[df$isfwd == -1,]);

	for( val in binVals)
	{
		#Get no of fwd prots that fall within bin and same for rev prots
		nFbin = nrow(df[df[,binCol]==val & df$isfwd !=-1 & df$decoyProb <= 0.5,]);
		nRbin = nrow(df[df[,binCol]==val & df$isfwd ==-1,]);

		#Calculate probabilities
		pgvnF = nFbin / nF;
		pgvnR = nRbin / nR;

		#Fill calculated probabilities in column
		df$pgvnFcol[df[,binCol]==val] = pgvnF;
		df$pgvnRcol[df[,binCol]==val] = pgvnR;
	}


	#Calculate adjusted probabilities and store in a vector
	df$adjProbsCol = (df$pgvnFcol*(1 - df$decoyProb) )/( df$pgvnFcol*(1 - df$decoyProb) + df$pgvnRcol*df$decoyProb ) ;

#print(df[df$adjProbsCol > 1,c("protid","isfwd","maxiniprob","decoyProb","pgvnFcol","pgvnRcol","adjProbsCol")]);

	return(df[,c(binCol,"pgvnFcol","pgvnRcol","adjProbsCol")]);
}

getAdjProbCol = function(df, binningCol, binThresholds)
{
	#Mark bins for the prots based on their binningCol vals
	#Store the results in adjDf because the new cols that 
	#will be added will be added again for the next iteration
	adjDf = markBins(df, binningCol, binThresholds);

	#Recreating the col name that was added in the markBins function
	binCol = paste(binningCol,"Bins",sep="");

	#Now calculate the adjusted True identification probabilty
	adjProbDf = calcAdjProb(adjDf, binCol);

	return(adjProbDf);
}


#Function takes df, does a random assignment of suppl values to the
#decoys and computes adjP for all proteins based on this assignment.
#The process is then repeated niter times and adjP values from all the
#iterations are collected into two matrices that are returned.

getAdjProbMatrix = function(df, sp, targetFdrs)
{
	niter=500;
	nrowsDf = nrow(df);

	rpkmBinThresholds = c(-8,-6,-4,-2,0,2,4,6,8);
	gpmBinThresholds = c(0,2,4,6,8,10,12,14);
	
	#Sort in increasing order of protlen
	#(Sorted dataframes are needed within the assignSupplValsToRev fn.
	#Doing it outside and passing the sorted dataframe saves time by
	#not having to do it every iteration)
	df = df[order(df$protlen),];
	sp = sp[order(sp$protlen),];
	
	#Split up the dataframe into separate fwd and reverse dataframes
	#(Faster to send only the rev dataframe for assigning Suppl vals
	fwd = df[df$isfwd ==1,];
	rev = df[df$isfwd !=1,];

	boundsList = getSupplValsSamplingBounds(rev$protlen,sp$protlen);	

	#Create empty matrices in which we can store the adjustedProbs
	rAdjProbMatrix = matrix(data=0, ncol=niter, nrow=nrowsDf);
	gAdjProbMatrix = matrix(data=0, ncol=niter, nrow=nrowsDf);
	combAdjProbMatrix = matrix(data=0, ncol=niter, nrow=nrowsDf);
	
	for(i in c(1:niter))
	{
		print(paste("Iteration no:",i));

		#Assign suppl vals to the rev prots by sampling from sp 
		rev = assignSupplValsToRev(rev,sp,boundsList);

		#Stick the rev vals back with fwd
		df = rbind(fwd,rev);

		print("Rev vals assigned");

		#Add log cols
		df$logrpkm = log2(df$rpkm);
		df$loggpm = log2(df$gpmNobs);

		#Get thresholds Df - adjusted by logRPKM
		rAdjPDf = getAdjProbCol(df,"logrpkm",rpkmBinThresholds);

		names(rAdjPDf) = c("rBin","pRgvnF","pRgvnR","rAdjP");

		rAdjPcol = rAdjPDf$rAdjP;
		pRgvnF = rAdjPDf$pRgvnF;
		pRgvnR = rAdjPDf$pRgvnR;
		
		#Get thresholds Df - adjusted by logGPMNobs
		gAdjPDf = getAdjProbCol(df,"loggpm",gpmBinThresholds);

		names(gAdjPDf) = c("gBin","pGgvnF","pGgvnR","gAdjP");

		gAdjPcol = gAdjPDf$gAdjP;
		pGgvnF = gAdjPDf$pGgvnF;
		pGgvnR = gAdjPDf$pGgvnR;

		#Use the conditional probabilities returned in adjPDf to
		#calculate the combined adjusted probability
		combAdjPcol = (pRgvnF*pGgvnF*(1 - df$decoyProb)) /( (pRgvnF*pGgvnF*(1 - df$decoyProb)) + (pRgvnR*pGgvnR*(df$decoyProb)) );

		#Store the adjThDfs returned into the appropriate lists
		rAdjProbMatrix[,i] = rAdjPcol;
		gAdjProbMatrix[,i] = gAdjPcol;
        combAdjProbMatrix[,i] = combAdjPcol;

		## Following lines Just for testing purposes; delete later
		# Store the adjProb cols into the df
		dfr = data.frame(df,rAdjPDf,gAdjPDf,combAdjPcol);

		#Write out an intermediate outfile
#		write.table(dfr,file=paste("iter.",i,".out.tsv",sep=""),sep="\t",quote=F,row.names=F);
	}

	#Because of the random assignment of rpkm and gpm values, in some
	#iterations we might end up with cases where no decoys are assigned
	#values from a certain bin. So for that bin, pGvnR =0. When this
	#happens to a protein that also has decoyProb =1, the adjP will be
	# (0/0) = NA. In these cases, the decoyProb = 1 will take precedence
	#since it is a more absolute value (not by random assignment). So
	#we set all NA values to be 0.
	rAdjProbMatrix[is.na(rAdjProbMatrix)] = 0;
	gAdjProbMatrix[is.na(gAdjProbMatrix)] = 0;
	combAdjProbMatrix[is.na(combAdjProbMatrix)] = 0;

	#Include protid and isfwd columns in the adjProbMatrix so that it
	#is easy to keep track of which value is for which protein
	rAdjProbMatrix = cbind(df[,c("protid","isfwd")],rAdjProbMatrix);
	gAdjProbMatrix = cbind(df[,c("protid","isfwd")],gAdjProbMatrix);
	combAdjProbMatrix = cbind(df[,c("protid","isfwd")],combAdjProbMatrix);

	allMatrixList = list(rAdjProbMatrix,gAdjProbMatrix,combAdjProbMatrix);

	return(allMatrixList);
}

calcFinalAdjP = function(df,colNames)
{
	#Calculate mean of the adjP values
	df$mean = apply(df[,colNames],1,mean);

	#Copy over as final cols
	df$final = df$mean;

	#Get avgDecoyDist

	#create a blank array to store avgDecoyDist
	avgDecoyDist = rep(0,nrow(df[df$isfwd==0,]));

	#Get the Sorted sum of adjP values, for all decoys
	#(Sorted sum - just sort the values and add all the first
	#rank values together, second rank values together etc.
	#The idea behind this is that, it isn't the actual decoy
	#value but just the overall dist of values that matters)
	for( col in colNames)
	{
		sortVals = sort(df[df$isfwd==0,col]);

		avgDecoyDist = avgDecoyDist + sortVals;
	}

	#Now that we have the sorted sum, divide by the no. of
	#iterations (also the no. of columns we added) to get the
	#average decoy distribution.
	avgDecoyDist = avgDecoyDist / length(colNames);

	#Sample values from the avgDecoyDist and fill in adjP final values
	#for the decoys
	df$final[df$isfwd ==0] = sample(avgDecoyDist,length(avgDecoyDist));

	return(df);
}


getFdrThresholds = function(df, sortCol, r, targetFdrs)
{

	#Just a hacky way to make decoyProb col get sorted ascending instead of descending
	if( sortCol== "decoyProb")
	{
		df$iDprob = 1 - df$decoyProb;
		sortCol = "iDprob";
	}

	#To start with, sort the df based on sortCol
	df = df[order(-df[,sortCol]),];

	#Get total no. of rows in df
	nrowsDf = nrow(df);

	#Defining empty data frames to rbind the result rows to later
	allFdrDf = data.frame(sortCol=numeric(), fdr=numeric(), nfwd=integer(),nrev=numeric());

	targetFdrDf = allFdrDf; #Just a quick way to define same empty df

	endIdx = 0;

	maxTargetFdr = tail(targetFdrs,n=1);

	while(endIdx < nrowsDf)
	{
		#Get sortVal at row next to endIdx
		#(If first iteration, get sortVal of row 1)
		sortVal = df[endIdx + 1, sortCol];

		#Now update endIdx value to the last row with sortCol
		#value greater than or equal to the value we are currently
		#looking at.
		endIdx = nrow(df[df[,sortCol] >= sortVal,]);


		#Copy the isfwd col into tmp arr
		tmp = df$isfwd[1:endIdx];

		#Get nfwd, nrev and calculate fdr
		fwd = sum(tmp == 1);
		rev = sum(tmp == 0);

		fdr = (r*rev) / fwd;

		if( fdr > maxTargetFdr & endIdx < nrowsDf)
		{
			#Find remaining no. of forward proteins in df
			nFwdRemain = sum(df$isfwd[(endIdx+1):nrowsDf] == 1);

			#Calculate the least possible value the Fdr can drop to after this point
			#(Assuming all the fwd proteins remaining in the df line up together 
			# after this point)
			minThFdr = (r*rev) / (fwd + nFwdRemain);

#			print(paste("f=",f,"r=",r,"fdr=",fdr,"nrowsDf=",nrowsDf,"endIdx=",endIdx,"nFwdRemain=",nFwdRemain,"minThFdr=",minThFdr));

			#If the minimum theoretical Fdr is greater than the maxTargetFdr,
			#we can break out of the loop because we can be sure that the fdr
			#will never drop lower than the last Fdr we want to find
			if( minThFdr > maxTargetFdr)
			{
				break;
			}
		}

		#Bind these values into allFdrDf
		allFdrDf = rbind(allFdrDf, c(sortVal, fdr, fwd, rev));
	}

	names(allFdrDf) = c(sortCol,"fdr","nfwd","nrev");

	for( fdrVal in targetFdrs)
	{
		#Extract the last row in allFdrDf that has fdr within fdrVal
		fdrValMaxRow = tail(allFdrDf[allFdrDf$fdr <= fdrVal,],n=1);

		if(nrow(fdrValMaxRow) ==0)
		{
			fdrValMaxRow[1,] = c(0,fdrVal,0,0);
		}

		targetFdrDf = rbind(targetFdrDf, fdrValMaxRow);
	} 

	#In targetFdrDf, replace the actual fdr values with the fdr 
	#threshold that they satisfy.
	targetFdrDf$fdr = targetFdrs;

	#Rename first column from 'sortCol' to the actual sort column name
	names(targetFdrDf) = names(allFdrDf);
	
	return(targetFdrDf);
}


makeImprovTable = function(df,r,targetFdrs)
{

 basicDf = getFdrThresholds(df,"maxhyperscore",r,targetFdrs);
 dprobDf = getFdrThresholds(df,"decoyProb",r,targetFdrs);
 rAdjDf = getFdrThresholds(df,"rfinal",r,targetFdrs);
 gAdjDf = getFdrThresholds(df,"gfinal",r,targetFdrs);
 combAdjDf = getFdrThresholds(df,"combfinal",r,targetFdrs);

 #Calculate percent improvement
 dImprov = (dprobDf$nfwd - basicDf$nfwd) / basicDf$nfwd;
 rImprov = (rAdjDf$nfwd - basicDf$nfwd) / basicDf$nfwd;
 gImprov = (gAdjDf$nfwd - basicDf$nfwd) / basicDf$nfwd;
 combImprov = (combAdjDf$nfwd - basicDf$nfwd) / basicDf$nfwd;

 #Put it all together as improvDf. Name cols appropriately.
 improvDf = data.frame(basicDf$fdr, basicDf$nfwd, rImprov, gImprov, combImprov, dImprov);
 names(improvDf)[1:2] = c("fdr","nfwd");

 #Also, include the rfactor used as a column. Will be helpful to manually
 #confirm the values later
 improvDf$rfactor = rep.int(r,nrow(improvDf));

 return(improvDf);
}


########### Main ############

wrapperFn = function(abacusOutfile,supplFile)
{
	targetFdrs = c(0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1);

	#Just extracting the file name from abacusOutfile
	abPrefix = unlist(strsplit(abacusOutfile,"/"))[5];

	#Read in abacusOutfile
	df = read.table(abacusOutfile, sep="\t",header=T,as.is=T);

	#AbacusOutfiles have headers in CAPS. Convert to lowercase
#	names(df) = tolower(names(df));

	#Read in supplFile
	sp = read.table(supplFile,sep="\t",header=T,as.is=T);

	#The abacusFile contains new.verbose output(containing all prots from
	#group, not just the representative protein). Reformat the df to
	#mark the repProts (not needed when using abacusout.mxvals.tsv as input)
#	print("Mark rep prots..");
#	df = markRepProts(df);

	#Read in decoyGroups file and use that information to mark the
	#decoys to use for decoyProb calculation
	print("Mark Dprob decoys..");
	df = markDprobDecoys(df,"~/vcap.rnaseq/input.files/decoyGroups.1.tsv");

	rd1 = calcRfactor(df,-1);

	#Calculate decoyProb for each protein (localised FDR)
	print("Calc decoy probs..");
	df = calcDecoyProb(df,rd1);
	
	#Add supplVals to the fwd prots in the Df.
	print("Add suppl fwd vals..");
	df = addSupplFwdVals(df,sp);

	#Add the decoyProb values to proteins in sp too (any prots not in df
	#get a decoyProb = 1)
	sp = addDprobToSuppl(sp,df);

	#We only kept the nonRep prots in order to add the DecoyProbs to 
	#the suppl dframe. For the rest of the processing and calculating
	#fdr thresholds we only need the rep Prots.
	df = df[df$isrep == 1,];
	
	#Now assign suppl vals to rev and get adjusted probabilities
	#This process will be repeated niter times and adjusted from
	#all the iterations will be returned in the form of a matrix.
	print("get thresholds Dfs list");
	allMatrixList = getAdjProbMatrix(df,sp,targetFdrs);

	#Separate out the two matrices from the lists
	rAdjMatrix = allMatrixList[[1]];
	gAdjMatrix = allMatrixList[[2]];
	combAdjMatrix = allMatrixList[[3]];

	#Write out the matrix of values to intermediate outfiles
	write.table(rAdjMatrix,file=paste(abPrefix,".r.tsv",sep=""),sep="\t",row.names=F,quote=F);
	write.table(gAdjMatrix,file=paste(abPrefix,".g.tsv",sep=""),sep="\t",row.names=F,quote=F);
	write.table(combAdjMatrix,file=paste(abPrefix,".comb.tsv",sep=""),sep="\t",row.names=F,quote=F);

	#Use tha adj values matrix to calculate a final adjusted value
	rFinal = calcFinalAdjP(rAdjMatrix,names(rAdjMatrix[,-c(1:2)]));
	gFinal = calcFinalAdjP(gAdjMatrix,names(gAdjMatrix[,-c(1:2)]));
	combFinal = calcFinalAdjP(combAdjMatrix,names(combAdjMatrix[,-c(1:2)]));
	
	#Add the final adjusted values columns to original df

	###Note:If I can keep track of the row order and make sure they are same,
	###I can directly assign the values without merging. That will make it a bit faster
	df = merge(df,rFinal[,c("protid","final")],by="protid");
	df$rfinal =df$final;
	df$final = NULL;

	df = merge(df,gFinal[,c("protid","final")],by="protid");
	df$gfinal = df$final;
	df$final = NULL;

	df = merge(df,combFinal[,c("protid","final")],by="protid");
	df$combfinal = df$final;
	df$final = NULL;
	
	print("Names df");
	print(names(df));

	#After calculating the adj probabilities, we don't need the -1 decoys anymore
	#so I can remove them. (If possible I should remove them before the calcFinalAdjP itself.
	#But I guess leaving them probably doesn't affect speed all that much(??) )
	df = df[df$isfwd != -1,];	

	#Calculate the rFactor using 0 decoys
	rd0 = calcRfactor(df,0);

	#Create a table of % improvements at various FDRs from df. (Compute fdr using the calculated r value)
	improvDf = makeImprovTable(df,rd0,targetFdrs);

	#Write out the df and improvDf to outfiles
	write.table(df,file=paste(abPrefix,".out.tsv",sep=""),sep="\t",row.names=F,quote=F);
	write.table(improvDf,file=paste(abPrefix,".improv.tsv",sep=""),sep="\t",row.names=F,quote=F);

}

##Main##
#abacusOutfiles = c("~/tpp.tmp/VCaP/fulldb/v40/v40.new.verbose.abacusout.tsv","~/tpp.tmp/VCaP/fulldb/v40/v40.new.verbose.abacusout.tsv","~/tpp.tmp/VCaP/fulldb/v40/v40.new.verbose.abacusout.tsv","~/tpp.tmp/VCaP/fulldb/v400/v400.new.verbose.abacusout.tsv","~/tpp.tmp/VCaP/fulldb/v40/v40.new.verbose.abacusout.tsv");
#abacusOutfiles = c("/tpp.data/avinashs/rna-seq/hek293/basic/samples2/h30.10/h30.10.s2.new.verbose.abacusout.mxvals.tsv");

abacusOutfiles = c("~/tpp.tmp/VCaP/fulldb/samples/v20.10/v20.10.new.verbose.abacusout.mxvals.tsv");

#suppl.file = "~/vcap.rnaseq/suppl.files/hek293/combined/hek293.comb.suppl.tsv";
suppl.file = "~/vcap.rnaseq/suppl.files/vcap.ens66.m1G/vcap.fwd.suppl.tsv";

outfile = "~/vcap.rnaseq/bayesPipelineResults/v400.1.pimprov.tsv";

for(abacusOut.file in abacusOutfiles)
{
	print(paste("Running fn for file",abacusOut.file));

	Rprof("Rprof.out");

	wrapperFn(abacusOut.file,suppl.file);
	
	Rprof(NULL);
}

print("All done!");

