# 2022_mortality

>   Note: In all notebooks, you must first change the paths for retrieving data from `data/` and for saving figure files.

*   The natural mortality rate (cohort mortality and actuarial mortality) is explored in `naturalmoralityrate.nb`.  
*   The starvation mortality rate, as well as changes to vital rate assumptions detailed in the Appendices, are explored in `starvationmortalityrate.nb`.  
*   The contemporary PPMR is explored in `predprayrateio.jl`, and exports a csv file that is then imported into the subsequent mathematica notebooks.  
*   The predation mortality rate is explored in `predationrateA.nb`, `predationrateB.nb`, and `predationrateC.nb`. Note that `predationrateA.nb` exports numerical data exploring different predation-prey mass ratios, which `predationrateB.nb` imports to build the figure. `predationrateC.nb` explores the generalized predator scenario.  
*   The harvest mortality rate is explored in `harvestmortalityrate.nb`.  