#Calculate various forms of genetic load.
#For revisions of Sanjak et al.
#Load is defined as deviation from optimum 
#fitness.  
#For the Gaussian fitness fxn, that optimum
#value is 1.
#We calculate the following loads:
#Total = 1 - fitness
#Fixed = 1 - (fitness calculated only using fixed variants)
#Segregating = 1 - (fitness calculated only using segregating variants)
#This is a custom plugin based on fwdpy.
#What we return is the mean value of each load over time.

from fwdpy.fwdpy cimport TemporalSampler,sampler_base,custom_sampler,singlepop_t,uint
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.memory cimport unique_ptr
from libc.math cimport sqrt,pow,exp

cdef struct load_values:
    unsigned generation
    double total
    double fixed
    double seg

#Functions for re-use below
cdef double gaussian_fitness(double P, double opt,double VS) nogil:
    return exp(-1.*pow(opt-P,2.)/(2.*VS))

cdef pair[double,double] haplotype_sums(const singlepop_t * pop, const size_t diploid) nogil:
    """
    Returns sum of s for each haplotype, only for seg variants
    """
    cdef pair[double,double] rv
    rv.first=0.
    rv.second=0.
    cdef unsigned twoN=2*pop.diploids.size()
    cdef size_t mut = 0
    for mut in range(pop.gametes[pop.diploids[diploid].first].smutations.size()):
        if pop.mcounts[pop.gametes[pop.diploids[diploid].first].smutations[mut]] < twoN:
            rv.first += pop.mutations[pop.gametes[pop.diploids[diploid].first].smutations[mut]].s
    for mut in range(pop.gametes[pop.diploids[diploid].second].smutations.size()):
        if pop.mcounts[pop.gametes[pop.diploids[diploid].second].smutations[mut]] < twoN:
            rv.second += pop.mutations[pop.gametes[pop.diploids[diploid].second].smutations[mut]].s
    return rv
        
#The following functions calculate the mean load for each of the 3 models.
#Notes:
#The fixed load is constant for everyone, so we just calculate it once.

cdef load_values additive_load(const singlepop_t * pop,const unsigned generation) nogil:
    cdef load_values rv
    rv.generation=generation
    rv.total=0.0
    rv.fixed=0.0
    rv.seg=0.0

    cdef size_t i = 0
    cdef unsigned twoN = 2*pop.diploids.size()
    #Fixed load
    cdef double sum_fixed_effects = 0.
    for i in range(pop.mutations.size()):
        if pop.mcounts[i]==twoN:
            sum_fixed_effects += pop.mutations[i].s
    #2*sum_fixed_effects b/c everyone is a homozygote
    #for a fixation
    cdef double fixed_w = gaussian_fitness(2.*sum_fixed_effects,0.,1.)
    rv.fixed=1.-fixed_w

    #Seg and total loads
    cdef pair[double,double] hapsums
    i=0
    cdef size_t mut=0
    for i in range(pop.diploids.size()):
        rv.total += (1.-pop.diploids[i].w)
        hapsums = haplotype_sums(pop,i)
        rv.seg += (1.-gaussian_fitness(hapsums.first+hapsums.second,0.,1.))
    rv.total /= <double>pop.diploids.size()
    rv.seg /= <double>pop.diploids.size()
    return rv;
