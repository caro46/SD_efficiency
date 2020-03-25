# R packages on cedar

```bash
module load  nixpkgs  gcc/7.3.0 r/3.6.1
#when doing module spider r/3.6.1 suggests a specific nixpkgs version but only one can be loaded ....
export R_LIBS=~/R/x86_64-pc-linux-gnu-library/3.6/

R
>install.packages("ggplot2")
>install.packages("dplyr")
>install.packages("qtl")
>install.packages("onemap")

# did it by steps to avoid cedar failing
```
# Linkage
Our crosses are usually outcrosses so need to use `Onemap` followed by `R/qtl`.
