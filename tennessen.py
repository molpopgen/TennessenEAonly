from __future__ import print_function   
import fwdpy as fp
import fwdpy.fitness as fpw
import fwdpy.qtrait as qt
import fwdpy.demography as dem
import numpy as np
import pandas as pd
import getopt
import sys
import math

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
    n.extend(dem.exponential_size_change(9300,51200,205)) #E5
    return np.array(n,dtype=np.uint32)

#        opts,args = getopt.getopt(sys.argv[1:],"m:l:r:o:d:s:",["model=","sampler=","cores=","batches=","seed=","usage"])
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
    print ("\t--model string (gbr) : Trait value model must be one of gbr, additive, or multi")
    print ("\t--sampler string (None) : Must be one of VA or stats")
    print ("\t--cores int>0 (64) : Number of populations to simulate simultaneously using different threads")
    print ("\t--batches int>0 (1) : Number of sets of simulations to do.  Total # sims will be (# cores)*(# batches)")
    print ("\t--seed int (0) : Random number seed")
    print ("\t--usage : print this info")
    
def setup_fitness(fitnessString):
    """
    Return the fitness model to use
    """
    if fitnessString == 'gbr':
        return fpw.SpopGBR()
    elif fitnessString == 'additive':
        return fpw.SpopAdditive()
    elif fitnessString == 'multi':
        return fpw.SpopMult()
    else:
        print("fitness model must be gbr, additive, or multi")
        usage()
        sys.exit(0)

def setup_sampler(samplerString,len):
    """
    Return the sampler to apply during the simulation
    """
    if samplerString == 'VA':
        return fp.VASampler(len)
    elif samplerString == "stats":
        #The second argument is the phenotypic optimum, which is 0 in these sims
        return fp.QtraitStatsSampler(len,0.0)
    else:
        print("sampler type must be VA or stats")
        usage()
        sys.exit(0)

def write_output(sampler,output,REPID):
    """
    Get the data from each sampler, add
    a "replicate ID number" to each element,
    then write to an HDF5 file.
    """
    df=[pd.DataFrame(i) for i in sampler.get()]
    for i in df:
        i['rep']=[REPID]*len(i.index)
        REPID+=1
    if isinstance(sampler,fp.VASampler):
        output.append('cumVA',pd.concat(df))
    elif isinstance(sampler,fp.QtraitStatsSampler):
        output.append('popstats',pd.concat(df))
    else:
        raise RuntimeError("uh oh: sampler type not recognized for output.  We shouldn't have gotten this far!")

    return REPID

def main():

    ##Parse options
    try:
        opts,args = getopt.getopt(sys.argv[1:],"m:l:r:o:d:s:",["model=","sampler=","cores=","batches=","seed=","usage"])
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

    #Now, we can get to work
    REPID = 0
    N=7310

    nregions=[] #Simulate no neutral variants
    sregions=[fp.ExpS(0,1,1,L,dominance)] #The dominance term will be ignored for the GBR model
    recregions=[fp.Region(0,1,1)] #Recombination uniform along region
    
    nosampler = fp.NothingSampler(ncores)
    rng = fp.GSLrng(seed)
    nlist=get_nlist()
    hdfout = pd.HDFStore(outfile,'w',complevel=6,complib='zlib')
    for i in range(nbatches):
        sampler = setup_sampler(samplerString,ncores)

        pops = fp.SpopVec(ncores,N)
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
                                                 tsample, #This will end up not doing anything...
                                                 sigE)
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
                                                 tsample, #This will end up not doing anything...
                                                 sigE)
        #If tsample is such that the last generation would not get processed,
        #then process it so that the final generation is included
        if float(len(nlist))%float(tsample) != 0.0:
            fp.apply_sampler(pops,sampler)
            
        REPID=write_output(sampler,hdfout,REPID)

    hdfout.close()
if __name__ == "__main__":
    main()
