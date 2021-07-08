---
title: "An R package for analyses of bioassays and probit graphs"
tags:
- R
- Dose-response
- Bioassays
- Probit analysis
- Exposure tests
date: "`r format(Sys.time(), '%d %B %Y')`"
output:
  md_document: default
  pdf_document:
    latex_engine: xelatex
authors:
- name: Piyal Karunarathne^[first author]
  orcid: 0000-0002-1934-145X
  affiliation: 1
- name: Nicolas Poquet^[co-first author]
  orcid: 0000-0003-3928-6803
  affiliation: 2
- name: Pascal Milesi^[co-last author]
  orcid: 0000-0001-8580-4291
  affiliation: 1
- name: Pierrick Labbé^[co-last author]
  orcid: 0000-0003-0806-1919
  affiliation: “3,4”
bibliography: paper.bib
affiliations:
- name: Plant Ecology and Evolution, Department of Ecology and Genetics, Evolutionary
    Biology Centre and SciLifeLab, Uppsala University, Uppsala, Sweden
  index: 1
- name: Institut Pasteur de Nouvelle-Calédonie, URE-Entomologie Médicale, Nouméa,
    New Caledonia
  index: 2
- name: Institut Universitaire de France, 1 Rue Descartes, 75231 Cedex 05, Paris.
  index: 3
- name: Institut des Sciences de l’Evolution de Montpellier (UMR 5554, CNRS-UM-IRD-EPHE),
    Université de Montpellier, Montpellier, 34095 Cedex 5, France
  index: 4
---

# Summary

Dose-response relationships (also known as exposure–response relationships) reflect the effects of a substance (most of the time a xenobiotic or a chemical) on organisms (populations, tissues or cells). Dose-response analyses are widely used in broad research areas, from medicine and physiology to vector control and pest management in agronomy. Further, reporting the response of organisms to stressors is an essential component of many public policies (e.g., health, environment).
An ideal example is monitoring of resistance to xenobiotics. Since the 1950s, xenobiotics (e.g., insecti-, pesti-, fungicides) have been widely used to control populations of vectors or pests. As a response, resistance mechanisms have been selected in targeted populations, undermining their efficiency. Establishing and comparing the resistance level of various populations to various xenobiotics is at the core of world health organization (WHO) recommendations in order to define / adjust vector control strategies. It is usually done by exposing batches of individuals (adults or larvae) to varying doses of the xenobiotic to assess their responses (mortality or knock-down effect). Despite the availability of statistical approaches for such analyses, there had been a lack of easily accessible analytical infrastructure for it (the traditionally-used software Probit ran in Basic and several labs kept an old computer for it). In 2013, we developed an R script with a robust statistical background to ease the dose-mortality relationship analysis. It has been used in many studies (e.g., @alout2016; pocquet2014; @badolo2019; @assogba2016; @yameogo2021; @epelboin2021; @perrier2021), and is now recommended as good practice by the ANSES (the French national agency for health and environment safety). In order to make it even more user friendly, we have now developed it into an R package called ‘BioRssay’ with more flexibility and improved presentation of results.

# Statement of need

‘BioRssay’ is a comprehensive compilation of scripts in R language[@core2020] designed to analyze dose-response relationships (or exposure-response: mortality, knock-down effect, etc.)  from bioassays of one or more strains, lines, populations (but also, cells etc.). This package provides a complete analytic workflow from data quality assessment to statistical analyses and data visualization. In the first steps, base-mortality in the controls (i.e., mortality linked to the experiment itself, not the exposure) is taken into account by adjusting the data following Abott’s correction [@abbott1925]. Data are then analyzed using a generalized linear model (probit-link function) to generate mortality-dose regressions (which take over-dispersion into account and allow for mortality of 0 or 1). Linearity of the log-dose response for each population is then tested using a chi-square test between model predictions and observed data (significant deviations from linearity may reflect mixed populations for example). By default, doses lethal for 25%, 50% and 95% of the populations (LD25, LD50 and LD95 respectively) are computed with their 95% confidence intervals (CI), following @johnson2013 approach, which allows taking the heterogeneity of the data into account [@finney1971]. Probit analysis. Cambridge: Cambridge University Press. 350 p.). Otherwise, the user has the option to specify any LD level. Likelihood ratio tests (LRT) are then implemented to test for statistical significance of the differences in response between different populations/strains, and when necessary, Holm-Bonferroni method [@holm1979]. is performed to control the family-wise error inherent to multiple testing. Finally, the resistance ratios for the required LD(s) (RR25, RR50 and RR95, by default), i.e., the LD(s) for a given population divided by the LD(s) of the population with the lowest one (usually the susceptible reference), are calculated according to @robertson1992, with 95% confidence intervals. Customizable plots of the probit-transformed regressions are also drawn (e.g.,with or without the desired confidence intervals).

# Citations

# Acknowledgements

We want to thank Jérôme Chopard, Haoues Alout, Mylène Weill and Nicole Pasteur for their valuable comments on earlier versions of the script and its outputs.

# References


