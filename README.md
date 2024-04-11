## On the dynamics of mortality and the ephemeral nature of mammalian megafauna
### Taran Rallings, Christopher P. Kempes, Justin D. Yeakel

Author Contact Information:   
- T Rallings (trallings-at-gmail-dot-com)  
- CP Kempes (ckempes-at-santafe-dot-edu)  
- JD Yeakel (jyeakel-at-ucmerced-dot-edu)  

Abstract - Energy flow through consumer-resource interactions is largely determined by body size. Allometric relationships govern the dynamics of populations by impacting rates of reproduction, as well as alternative sources of mortality, which have differential impacts on smaller to larger organisms. Here we derive and investigate the timescales associated with four alternative sources of mortality for terrestrial mammals: mortality from starvation, mortality associated with aging, mortality from consumption by predators, and mortality introduced by anthropogenic subsidized harvest. The incorporation of these allometric relationships into a minimal consumer-resource model illuminates central constraints that may contribute to the structure of mammalian communities. Our framework reveals that while starvation largely impacts smaller-bodied species, the allometry of senescence is expected to be more difficult to observe. In contrast, external predation and subsidized harvest have greater impacts on the populations of larger-bodied species. Moreover, the inclusion of predation mortality reveals mass thresholds for mammalian herbivores, where dynamic instabilities may limit the feasibility of megafaunal populations. We show how these thresholds vary with alternative predator-prey mass relationships, which are not well understood within terrestrial systems. Finally, we use our framework to predict the harvest pressure required to induce mass-specific extinctions, which closely align with previous estimates of anthropogenic megafaunal exploitation in both paleontological and historical contexts. Together our results underscore the tenuous nature of megafaunal populations, and how different sources of mortality may contribute to their ephemeral nature over evolutionary time.


## Repository
Included here are the raw data and julia + mathematica notebooks required to generate the primary and supplementary figures in the manuscript.

## Data Files (/data/)
*   `densitydata.csv`: Mammalian density data by body size from Damuth et al. (1981). Columns - C1: Body mass (grams); C2: Density (inds/m^2)  
*   `carnivoredensities_trimmed.csv`: Mammalian carnivore density data by body size from Carbon & Gittleman (2002). Columns - C1: Body mass (grams); C2: Density (inds/m^2)  
*   `data_hayward_all.csv`: Prey preferences and dietary data for large-bodied carnivores collated from Hayward, 2006; Hayward et al., 2006a,b; Hayward and Kerley, 2008, 2005; Hayward et al., 2006c. Columns - C1: Predator; C2: Predbodymasskg; C3: Prey; C4: JacobsIndex; C5: JI_SE; C6: PercentOfKills; C7: SE; C8: Preybodymasskg34adultfemalemass; C9: Preybodymasskg; C10: HerdSize; C11: HabitatDensity; C12: ThreatToPredator. See specific references for details
*   `ppmr_fit_table_revreps.csv`: Predator-prey mass relationship calculated in the Julia file `predpreyratio.jl`. These are the intercepts and slopes that determine the expected predator mass given a prey mass, and are imported and used in the mathematica scripts. Columns - C1: Best Fit; C2: 95% confidence interval low; C3: 95% confidence interval high. Rows - R1: Intercept; R2: Slope.  
*   `predpreymass_tablereps.csv`: Bootstrapped values of predator and prey masses calculated in the Julia file `predpreyratio.jl`. These values are used to calculate the above predator-prey mass relationship. Columns - C1: prey mass (kg); C2: pred mass (kg)

## Scripts (/)
>   Note: In all scripts, you must first insert your own file path for retrieving data from `/data/` and for saving figure files.

### Julia (v. 1.8.2)
*   `predpreyratio.jl`: This Julia script imports the averaged data collated from Hayward, 2006; Hayward et al., 2006a,b; Hayward and Kerley, 2008, 2005; Hayward et al., 2006c, (`data_hayward_all.csv`) describing the expected prey size for different large-bodied mammalian carnivores. The script then resamples predator body masses for each prey to obtain the expected predator mass given a particular prey body size. The fitted relationship (and 5/95% confidence intervals) are exported as `ppmr_fit_table_revreps.csv`, and the raw bootstrapped samples are exported as `predpreymass_tablereps.csv` (see above). These arrays are imported into the various mathematica scripts below for analysis and figure creation. 

### Mathematica (v. 14.0)

*   `carryingcapacity_stability.nb`: Assessment of how changes to resource growth rates and carrying capacities alter the range of consumer mass-density relationships. This creates Figure 1 (main text). In addition, this script contains the stability analysis detailed in Supplemental Figure C1 (determinant of the Jacobian and real values of associated Jacobian eigenvalues for consumer-resource system as a function of consumer mass) and Supplemental Figure C5 (the effect of changes to the carrying capacity on the consumer mass threshold).  
*   `naturalmoralityrate.nb`: Details expression of the natural mortality rate (cohort mortality and actuarial mortality). This creates Figure 2 (main text).  
*   `starvationmortalityrate.nb`: Details expression of the the starvation mortality rate. This creates Figure 3 (main text).  
*   `predationrate_specialist.nb`: Details expression of the mortality rate due to a specialist predator featured primarily in the main text, and creates Figure 4 (main text). This script imports `ppmr_fit_table_revreps.csv` and `predpreymass_tablereps.csv` to build panel B in Figure 4, and establishes the expected predator body size given a prey body size used in the analysis. This script also creates Supplemental Figure C3 (mass ranges corresponding to feasible megatrophic interactions).
*   `predationrate_generalist.nb`:  Details expression of the mortality rate due to a generalist predator, which is discussed in the main text, with findings primarily shown in the Supplement. Creates Supplemental Figure C2 (the effect of changing the predation intensity on the single herbivore consumer population) and C4 (the effects of lower predation intensity threshold herbivore mass and threshold predator mass across variable PPMRs) in the supplement. 
*   `harvestmortalityrate.nb`: Details expression of the harvest pressure. This creates Figure 5 (main text).   


### References
*   Damuth, J. 1981. Population density and body size in mammals. Nature 290:699–700.  
*   Carbone, C., and J. L. Gittleman. 2002. A common rule for the scaling of carnivore density. Science 295:2273–2276.
*   Hayward, M. 2006. Prey preferences of the spotted hyaena (*Crocuta crocuta*) and degree of dietary overlap with the lion (*Panthera leo*). J. Zoology 270:606–614.  
*   Hayward, M., P. Henschel, J. O'Brien, M. Hofmeyr, G. Balme, and G. I. Kerley. 2006a. Prey preferences of the leopard (*Panthera pardus*). Journal of Zoology 270:298–313.  
*   Hayward, M., M. Hofmeyr, J. O'brien, and G. I. Kerley. 2006b. Prey preferences of the cheetah (*Acinonyx jubatus*)(felidae: Carnivora): morphological limitations or the need to capture rapidly consumable prey before kleptoparasites arrive? Journal of Zoology 270:615–627.  
*   Hayward, M. W., and G. Kerley. 2008. Prey preferences and dietary overlap amongst Africa's large predators. S. African J. Wild. Res. 38:93–108.  
*   Hayward, M. W., and G. I. Kerley. 2005. Prey preferences of the lion (*Panthera leo*). Journal of zoology 267:309–322.  
*   Hayward, M. W., J. O'Brien, M. Hofmeyr, and G. I. Kerley. 2006c. Prey preferences of the african wild dog *Lycaon pictus* (canidae: Carnivora): ecological requirements for conservation. Journal of Mammalogy 87:1122–1131  