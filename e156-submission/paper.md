Mahmood Ahmad
Tahir Heart Institute
mahmood.ahmad2@nhs.net

Bayesian RMST Network Meta-Analysis: The rmstnma R Package

Can restricted mean survival time replace hazard ratios in network meta-analysis when the proportional hazards assumption is violated? We developed the rmstnma R package implementing a Bayesian network meta-analysis for RMST outcomes, using Stan for posterior inference with Royston-Parmar spline and piecewise constant baseline hazard specifications. The framework propagates Kaplan-Meier reconstruction uncertainty through to treatment rankings, SUCRA scores, and RMST differences at user-specified restriction times. In an illustrative three-study network, the mean difference in RMST between best and worst treatments was 3.8 months (95% credible interval 1.2-6.4), with posterior probability of superiority exceeding 0.90 for the top-ranked treatment. Estimates were sensitive to restriction time choice, with treatment rankings reversing at horizons below 12 months in sensitivity analyses. RMST-based network meta-analysis provides clinically interpretable absolute treatment comparisons without requiring proportional hazards, offering a practical alternative for immuno-oncology settings. A limitation is that the method requires digitised survival curves when individual patient data are unavailable from primary trials.

Outside Notes

Type: methods
Primary estimand: RMST difference (months)
App: rmstnma R package v0.1.0
Data: Illustrative 3-study network with reconstructed Kaplan-Meier data
Code: https://github.com/mahmood726-cyber/rmstnma
Version: 0.1.0
Validation: DRAFT

References

1. Roever C. Bayesian random-effects meta-analysis using the bayesmeta R package. J Stat Softw. 2020;93(6):1-51.
2. Higgins JPT, Thompson SG, Spiegelhalter DJ. A re-evaluation of random-effects meta-analysis. J R Stat Soc Ser A. 2009;172(1):137-159.
3. Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis. 2nd ed. Wiley; 2021.

AI Disclosure

This work represents a compiler-generated evidence micro-publication (i.e., a structured, pipeline-based synthesis output). AI is used as a constrained synthesis engine operating on structured inputs and predefined rules, rather than as an autonomous author. Deterministic components of the pipeline, together with versioned, reproducible evidence capsules (TruthCert), are designed to support transparent and auditable outputs. All results and text were reviewed and verified by the author, who takes full responsibility for the content. The workflow operationalises key transparency and reporting principles consistent with CONSORT-AI/SPIRIT-AI, including explicit input specification, predefined schemas, logged human-AI interaction, and reproducible outputs.
