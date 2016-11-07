from __future__ import print_function   
import fwdpy as fp
import fwdpy.qtrait as qt
import fwdpy.demography as dem
import numpy as np
import pandas as pd
import getopt
import sys
import math
import datetime
import gc
#This is custom plugin:
import pyximport
pyximport.install()
import loads

def get_nlist():
    """
    Generates a numpy array of the canges in N over time

    There are 5 epochs, with t=0 being the present.
    
    E1: Ne= 7,310  from t=start(8N?) to t = - 5920 generation (Ancestral sizes until 5920 generations ago)
    
    E2: Ne =14,474 from t = -5920 to t = -2040 (Ancient growth at 5920 g ago)
    
    E3: Ne =1,861 from t= -2040 to t= -920 (OOA, first bottle neck 2040 g ago)
    
    E4: Ne = 1,032 to Ne = 9,300 during t = -920 to t = -205 ( second bottle neck and onset of 715 g of exponential growth at rate 0.31% per gen )  
    
    E5: Ne = 9,300 to Ne = 512,000 during t = -205 to t = -0 ( 205 g of exponential growth at rate 1.95% per gen )  
    """
    n=[7310]*(10*7310) #E1: evolve ancestral size to mutation/selection/drift equilibrium
    n.extend([14474]*(5920-2040)) #E2
    n.extend([1861]*(2040-920)) #E3
    n.extend(dem.exponential_size_change(1032,9300,920-205)) #E4
    n.extend(dem.exponential_size_change(9300,512000,205)) #E5
    return np.array(n,dtype=np.uint32)

def get_epoch_lengths():
    """
    Gets the epoch lengths of the demographic model
    """
    rv=[]
    rv.append(10*7310)
    rv.append(5920-2040)
    rv.append(2040-920)
    rv.append(920) #We merge the last two into a single "growth epoch" so that we only have 1 really big regression to do
    return rv

def usage():
    print ("Usage:")
    print ("python tennessen.py [options]")
    print ("Default option values shown in parentheses")
    print ("Defaults of None must be provide by user")
    print ("options are (defaults shown in parenthesis")
    print ("\t-m float>0 (1.25e-4) : a float specifying mutation rate to causal variants (per gamete, per generation)")
    print ("\t-l float>0 (None) : mean effect size of a causal variant")
    print ("\t-r float>0 (1.25e-3) : crossover rate (per diploid, per generation)")
    print ("\t-o string (None) : output file name, which will be in HDF5 format")
    print ("\t-d float (1.0) :  Dominance of causal mutations.  Only relevant to additive and multiplicative trait value models")
    print ("\t-s float>0 (0.075) : Standard deviation in environmental noise added to trait values.")
    print ("\t-t int>0 (50) : Apply the sampler every t generations")
    print ("\t--model string (gbr) : Trait value model must be one of gbr, additive, or multi")
    print ("\t--sampler string (None) : Must be one of VA, stats, or load")
    print ("\t--bigstub string (None) : Base file name for genotype matrices for N=512,000.  Required for VA sampler and ignored otherwise")
    print ("\t--cores int>0 (64) : Number of populations to simulate simultaneously using different threads")
    print ("\t--batches int>0 (1) : Number of sets of simulations to do.  Total # sims will be (# cores)*(# batches)")
    print ("\t--seed int (0) : Random number seed")
    print ("\t--usage : print this info")
    
def setup_fitness(fitnessString):
    """
    Return the fitness model to use
    """
    if fitnessString == 'gbr':
        return qt.SpopGBRTrait()
    elif fitnessString == 'additive':
        return qt.SpopAdditiveTrait()
    elif fitnessString == 'multi':
        return qt.SpopMultTrait()
    else:
        print("fitness model must be gbr, additive, or multi")
        usage()
        sys.exit(0)

def setup_sampler(samplerString,fitness,len):
    """
    Return the sampler to apply during the simulation
    """
    if samplerString == 'VA':
        return fp.VASampler(len)
    elif samplerString == "stats":
        #The second argument is the phenotypic optimum, which is 0 in these sims
        return fp.QtraitStatsSampler(len,0.0)
    elif samplerString == "load":
        if isinstance(fitness,qt.SpopGBRTrait):
            return loads.gbrLoad(len)
        elif isinstance(fitness,qt.SpopMultTrait):
            return loads.multiplicativeLoad(len)
        elif isinstance(fitness,qt.SpopAdditiveTrait):
            return loads.additiveLoad(len)
    else:
        print("sampler type must be VA, stats, or load")
        usage()
        sys.exit(0)


def write_output(sampler,output,REPID):
    """
    Get the data from each sampler, add
    a "replicate ID number" to each element,
    then write to an HDF5 file.
    """
    df=[pd.DataFrame(i) for i in sampler.get()]
    if isinstance(sampler,fp.VASampler):
        output.append('cumVA',pd.concat(df))
    elif isinstance(sampler,fp.QtraitStatsSampler):
        output.append('popstats',pd.concat(df))
    elif isinstance(sampler,loads.gbrLoad):
        print ("gbr")
    elif isinstance(sampler,loads.additiveLoad):
        print ("additive")
    elif isinstance(sampler,loads.multiplicativeLoad):
        print ("multiplicative")
    else:
        raise RuntimeError("uh oh: sampler type not recognized for output.  We shouldn't have gotten this far!")

    return REPID

#def write_files(args):
#def write_files(pops,REPID):
    #file_sampler=fp.GenoMatrixSampler(len(args[0]),True,True)
#    file_sampler=fp.GenoMatrixSampler(len(pops),True,True)
#    fp.apply_sampler(pops,file_sampler)
#    file_sampler.tofile("foobig",REPID)
#    file_sampler.force_clear()
#    del file_sampler
    
def main():
    ##Parse options
    try:
        opts,args = getopt.getopt(sys.argv[1:],"m:l:r:o:d:s:t:",["model=","sampler=","cores=","batches=","seed=","usage","bigstub="])
    except getopt.GetoptError as err:
        print(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    ##Establish default parameters
    sampler = None
    samplerString = None
    outfile = None
    fitness = None
    fitnessString = 'gbr'
    L = None #This is "lambda" = mean effect size
    ncores = 64 #default number of populations to simulate at once
    nbatches = 1 #default number of batches of ncores populations to simulate
    tsample = 50 #Default time interval to apply sampler
    seed = 0
    #These are the params from TFL2013,JSS et al
    mutrate = 1.25e-4
    recrate = 1.25e-3
    sigE = 0.075
    bigstub = None    
    #For "site-based" models with dominance
    dominance = 1.0
    
    for o,a in opts:
        if o == '-m':
            mutrate=float(a)
            if mutrate <= 0.0:
                print("mutation rate must be positive")
                usage()
                sys.exit(0)
        elif o == '-l':
            L = float(a)
        elif o == '-r':
            recrate=float(a)
            if recrate <= 0.0:
                print("recombination rate must be positive")
                usage()
                sys.exit(0)
        elif o == '-o':
            outfile = a
        elif o == '-d':
            dominance = float(a)
        elif o == '-s':
            sigE=float(a)
            if sigE < 0.0:
                print("sigmaE must be >= 0.0")
                usage()
                sys.exit(0)
        elif o == '-t':
            tsample = int(a)
            if tsample < 1:
                print("sampling interval must be > 0")
                usage()
                sys.exit(0)
        elif o == '--model':
            fitnessString = a
        elif o == '--sampler':
            samplerString = a
        elif o == '--cores':
            ncores = int(a)
            if ncores<1:
                print("number of cores to use must be >= 1")
                usage()
                sys.exit(0)
        elif o == '--batches':
            nbatches=int(a)
            if nbatches<1:
                print("number of batches to simulations must be >= 1")
                usage()
                sys.exit(0)
        elif o == '--seed':
            seed = int(a)
        elif o == '--bigstub':
            bigstub=a
        elif o == '--usage':
            usage()
            sys.exit(0)

    fitness = setup_fitness(fitnessString)        

    if L is None:
        print("mean effect size must be defined")
        usage()
        sys.exit(0)
    if outfile is None:
        print("output file undefined")
        usage()
        sys.exit(0)
    if samplerString=="VA" and bigstub is None:
        print("must define --bigstub when using VA sampler")
        usage()
        sys.exit(0)
    #Now, we can get to work
    REPID = 0
    N=7310

    nregions=[] #Simulate no neutral variants
    sregions=[fp.ExpS(0,1,1,L,dominance)] #The dominance term will be ignored for the GBR model
    recregions=[fp.Region(0,1,1)] #Recombination uniform along region
    
    nosampler = fp.NothingSampler(ncores)
    rng = fp.GSLrng(seed)
    nlist=get_nlist()
    epochs=get_epoch_lengths()
    hdfout = pd.HDFStore(outfile,'w',complevel=6,complib='zlib')
    for i in range(nbatches):
        #print ("starting batch",i,"of",nbatches," at",datetime.datetime.now().time().isoformat())
        gc.collect()
        pops = fp.SpopVec(ncores,N)
        if samplerString != 'VA':
            #We will evolve the first 8N generations w/o sampling anything
            qt.evolve_regions_qtrait_sampler_fitness(rng,
                                                     pops,
                                                     nosampler,
                                                     fitness,
                                                     nlist[:8*N+1],
                                                     0.0, #No neutral mutations
                                                     mutrate,
                                                     recrate,
                                                     nregions,
                                                     sregions,
                                                     recregions,
                                                     N, #This will end up not doing anything...
                                                     sigE)
            #print("evolved batch",i,"to equilibrium at",datetime.datetime.now().time().isoformat())
            sampler = setup_sampler(samplerString,fitness,ncores)
            #Now, evolve the rest of the way and sample...
            qt.evolve_regions_qtrait_sampler_fitness(rng,
                                                     pops,
                                                     sampler,
                                                     fitness,
                                                     nlist[(8*N+1):],
                                                     0.0, #No neutral mutations
                                                     mutrate,
                                                     recrate,
                                                     nregions,
                                                     sregions,
                                                     recregions,
                                                     tsample, 
                                                     sigE)
            #If tsample is such that the last generation would not get processed,
            #then process it so that the final generation is included
            if float(len(nlist))%float(tsample) != 0.0:
                fp.apply_sampler(pops,sampler)
                #print ("Applied final sampler at",datetime.datetime.now().time().isoformat())
                REPID=write_output(sampler,hdfout,REPID)
        else:
            #For the VA sampler, we simulate each epoch separately.
            #After 10N generations to equilibrium, we apply a VA
            #sampler at 1st & last generation of each epoch.
            #We treat both periods of growth as 1 epoch.
            se=0
            EPOCH=0
            for e in epochs:
                start=se
                end=se+e-1
                #print(e,se,len(nlist),nlist[start],nlist[end],len(nlist[start:end+1]))
                if EPOCH > 0:
                    #Evolve pop a single generation at beginning of this epoch
                    qt.evolve_regions_qtrait_sampler_fitness(rng,
                                                             pops,
                                                             nosampler,
                                                             fitness,
                                                             nlist[start:start+1],
                                                             0.0, #No neutral mutations
                                                             mutrate,
                                                             recrate,
                                                             nregions,
                                                             sregions,
                                                             recregions,
                                                             N, 
                                                             sigE)
                    sampler=fp.VASampler(len(pops))
                    fp.apply_sampler(pops,sampler)
                    dummy=write_output(sampler,hdfout,REPID)
                #Do the evolution for rest of epoch
                qt.evolve_regions_qtrait_sampler_fitness(rng,
                                                         pops,
                                                         nosampler,
                                                         fitness,
                                                         nlist[start+1:end+1],
                                                         0.0, #No neutral mutations
                                                         mutrate,
                                                         recrate,
                                                         nregions,
                                                         sregions,
                                                         recregions,
                                                         N, 
                                                         sigE)
                if EPOCH > 0:
                    #print ("ending epoch",EPOCH," of",len(epochs))
                    #We just don't have the
                    #RAM to do the regression in memory
                    #for this pop size
                    if EPOCH < len(epochs)-1:
                        sampler=fp.VASampler(len(pops))
                        dummy=write_output(sampler,hdfout,REPID)
                    else:
                        #We gotta do something else,
                        #which means write the damn things to file!
                        #print ("sampling files in epoch ",EPOCH)
                        file_sampler=fp.GenoMatrixSampler(len(pops),True,True)
                        fp.apply_sampler(pops,file_sampler)
                        file_sampler.tofile("foobig",REPID)
                        file_sampler.force_clear()
                        del file_sampler
                se+=e
                EPOCH+=1
        REPID+=len(pops)
        #print ("finished evolving batch",i,"at",datetime.datetime.now().time().isoformat())
    hdfout.close()
if __name__ == "__main__":
    main()
