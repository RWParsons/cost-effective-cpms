# "Cost Effective Clinical Prediction Models" (cost-effective-cpms)

This repository contains the code used to evaluate clinical prediction models in terms of Net Monetary Benefit under different cutpoint selection methods. This included the "value-optimizing" cutpoint method introduced in our publication here (<https://doi.org/10.1093/jamia/ocad042>).

> Integrating economic considerations into cutpoint selection may help align clinical decision support towards value-based healthcare.<br>
Rex Parsons, Robin Blythe, Susanna Cramb, Steven McPhail <br>
_Journal of the American Medical Informatics Association_ 2023

The repository already includes the outputs of the code but to repeat the analyses, the `./analyses.Rmd` file can be run.

In the manuscript, a shiny app is presented that can visualise selected cutpoints, including the proposed value-optimizing cutpoint (when NMB values are provided by the user). This is available here (<https://aushsi.shinyapps.io/cost-effective-cpms-app/>) and the code used to develop that app is available here (<https://github.com/RWParsons/cost-effective-cpms-app>).

Following on from this study, we developed `{predictNMB}`: an R package that encapsulates the same simulation process used in this study. It can be used to evaluate hypothetical clinical prediction models when they are used to guide treatment decisions, and evaluation is performed in terms of Net Monetary Benefit (NMB).

More details about `{predictNMB}` can be found on it's website (<https://rwparsons.github.io/predictNMB/>). <a href='https://rwparsons.github.io/predictNMB/'><img src='./hex.png' align="right" height="150" /></a>
