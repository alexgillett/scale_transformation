# scale_transformation
Demonstration of equations transforming effect size estimates between the logit and the liability scale.
Example R code to accompany the paper: 
Transforming summary statistics from logistic regression to the liability scale: application to genetic and environmental risk scores.

An important part of stratified medicine is the development of models of disease risk that are easily updated using model parameter estimates from the existing literature. That is, developing models of risk that can be parameterised without raw data. This will require parameter estimates from different studies to be combined into one model, and therefore a common scale for risk is required. Two important scales in genetic association studies are the logit and the liability scale. We therefore have derived approximations to translate univariate effect size estimates for risk variables from the logit to the liability scale, and vice versa. This allows researchers to select a common scale, and transform any parameter estimates to this scale if required.

In particular, we focus on transforming effect size estimates from the logit to the liability scale. This is because we wish to include a polygenic risk score (PRS) in the developed risk model, and a useful and commonly calculated summary measure for the PRS; the proportion of variability in liability to disease attributable to the PRS, is presented on the liability scale.

In this example script we demonstrate how to transform effect size estimates from the logit scale to the liability scale. This is done using schizophrenia as an example. We find effect size estimates for 5 environmental risk factors on the liability scale using odds ratios reported in the literature. We create an environmental risk score (ERS) using these transformed effect size estimates and then estimate the risk of schizophrenia using this ERS and an existing schizophrenia PRS.

