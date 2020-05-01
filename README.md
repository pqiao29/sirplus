# sirplus

sirplus is a package for the modeling of COVID-19 spread using stochastic individual compartment models (ICMs). The model implemented extends the classical Susceptible-Infectious-Recovered (SIR) model by adding compartments for Exposed, Quarantined, Hospitalised and Case Fatality individuals. The package provides a simple interface for creating experiments to demonstrate how factors like social-distancing and epidemiological parameters will change the curve.

## Installation and vignette

To view the vignette install the package then: 

```
devtools::install(build_vignettes=TRUE)
library(sirplus)
browseVignettes("sirplus")  
```

To view the vignette of the package, visit [the package website](https://pqiao29.github.io/sirplus/)
