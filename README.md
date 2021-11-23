
# kerneval

<img src="man/figures/kerneval_hex_sticker.png" width="215px" align="right">

The goal of the *kerneval* R package is to estimate probability densities from data that are affected by observation bias, implementing statistical theory from Jones (1991), Barmi and Simonoff (2000), and Borrajo et al. (2017). Many types of observation bias can be treated in the same way mathematically. As one example, length bias was first described by David Cox (1969) using data from textile manufacturing, but the same selection bias structure occurs in ecology with the chance of identifying an animal being directly proportional to group size. I originally developed *kerneval* to analyze fossil data sampled patchily over environments (Antell et al. 2021), but the functions are non-specific to paleontology.

If you are new to kernel density estimation of probability densities, the introductory sections of the vignette (`vignettes/kerneval_vignette`) explain the basic theory and use cases. The rest of the vignette reviews the relevant principles of bias correction and walks through several code examples. The `kerneval_technical_details` document elaborates on the methods to approximate integrals, apply boundary reflection, and calculate mean integrated squared error.

## Installation

You can install kerneval in R with:

``` r
devtools::install_github('GwenAntell/kerneval')
```

## References

Antell, G.S., Fenton, I.F., Valdes, P., and Saupe, E.E. (2021). "Thermal niches of planktonic foraminifera are static throughout glacial-interglacial climate change." Proceedings of the National Academy of Sciences. https://doi.org/10.1073/pnas.2017105118

Barmi, H.E., and Simonoff, J.S. (2000). "Transformation-Based Density Estimation for Weighted Distributions." Journal of Nonparametric Statistics 12 (6): 861–78.

Borrajo, M.I., González-Manteiga, W., and Martínez-Miranda, M.D. (2017). "Bandwidth Selection for Kernel Density Estimation with Length-Biased Data." Journal of Nonparametric Statistics 29 (3): 636–68.

Cox, D. (1969). "Some Sampling Problems in Technology." In: New Developments in Survey Sampling, edited by NL Johnson and H Smith, 506–27. New York: Wiley.

Jones, M.C. (1991). "Kernel Density Estimation for Length Biased Data." Biometrika 78 (3): 511–19.
