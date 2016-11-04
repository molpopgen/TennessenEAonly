#!/usr/bin/env Rscript
library(data.table)
library(biglm)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)

fillvpv1matrix <- function(mm,twoN)
    {
        rv = data.frame((1:twoN)/twoN,NA,NA)
        colnames(rv) <- c("V1","V2","V3")
        for( i in 1:nrow(mm) )
            {
                P = mm[i,1]*twoN
                rv[P,2]=mm[i,2]
                rv[P,3]=mm[i,3]
            }
        ##Courtesy of http://stackoverflow.com/questions/23340150/using-dplyr-window-functions-to-make-trailing-values
        rv %>%
            mutate(dummy = cumsum(0 + !is.na(V2))) %>%
                group_by(dummy,add=TRUE) %>%
                    mutate(filled = nth(V2,1)) %>%
                        mutate(filled2 = nth(V3,1)) %>%
                        ungroup() %>%
                            select(-V2,-V3,-dummy) -> m2
        colnames(m2) <- c("p","rsq","adj.rsq")
        return(m2)
    }

dobiglm<-function(data)
{
    #we have to get the allele counts as we read the data to avoid reading it back in later
    #mutc = colSums(data[,2:ncol(data)])
    #set the count of the dominance component to the count of the corresponding additive component
    #mutc[seq(2,length(mutc),2)] = mutc[seq(1,length(mutc)-1,2)]

    START=1
    STOP=5000
    subset=rep(0,popsize)

    subset[START:STOP]=1
    FORMULA <- paste(colnames(data)[1], "~", paste(colnames(data)[-1], collapse=" + "))

    fit=biglm(as.formula(FORMULA),data=data[which(subset==1),])
    subset[START:STOP]=0
    START=STOP+1
    STOP=min(popsize,START+5000)
    subset[START:STOP]=1
    while( START < popsize )
    {
        print(paste(START,STOP))
        update(fit,moredata=data[which(subset==1),])
        subset[START:STOP]=0
        START=STOP+1
        STOP=min(popsize,START+5000)
        subset[START:STOP]=1
    }
    #RV =list(SS=c(BIGGIE$qr$D*BIGGIE$qr$thetab^2,BIGGIE$qr$ss)[-1],ac=alleleCounts)
    RV=summary(aov(fit,data=data))
    return(RV)
}

data=fread(paste("gunzip -c ",args[1]))
popsize=nrow(data)

data.aov.s=dobiglm(data)
ROWS=rownames(data.aov.s[[1]])
print(ROWS)
alleleCounts = as.integer(colSums(data[,sapply(as.array(ROWS[1:(length(ROWS)-1)]),function(x) gsub(" ","",x),USE.NAMES=FALSE)]))

##Initalize the matrix to return.
##The columns will be: mutant allele frequency, r^2, adj. r^2
p = array( dim = length(unique(alleleCounts)) )
rsq = array( dim = length(unique(alleleCounts)) )
adj.rsq = array( dim = length(unique(alleleCounts)) )

sum.sq = data.aov.s[[1]]$'Sum Sq'
ttl.sum.sq = sum(sum.sq)
DF = data.aov.s[[1]]$'Df'
n = length(sum.sq)
##Populate the matrix, starting with the rares
IDX=1
for( ac in unique(sort(alleleCounts)) )
    {
        ac.sum.sq = sum.sq[which(alleleCounts == ac)]
        p[IDX] = ac/twoN
        ## r^2 due just to mutation at this freq. bin
        rsq[IDX] = sum(ac.sum.sq)/ttl.sum.sq
        if(is.na(rsq[IDX]))
            {
                stop("rsq == NA found")
            }
        ## adj. r^2 due just to mutations at this freq. bin
        adj.rsq[IDX] =  1 - ( (sum.sq[n] + sum(sum.sq[which(alleleCounts != ac)]))/ttl.sum.sq )*(sum(DF)/(DF[n] + length(which(alleleCounts != ac))))
        if(is.na(adj.rsq[IDX]))
            {
                stop("adj.rsq == NA found")
            }
        IDX=IDX+1
    }
rsq = cumsum(rsq)
adj.rsq = cumsum(adj.rsq)
final=fillvpv1matrix(data.frame(p,rsq,adj.rsq),2*popsize)
