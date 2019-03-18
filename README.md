# scale_transformation
Demonstration of equations transforming effect size estimates between the logit and the liability scale.
Example R code to accompany the paper: 
'Transforming summary statistics from logistic regression to the liability scale: application to genetic and environmental risk scores' [DOI: 10.1159/000495697].

An important part of stratified medicine is the development of models of disease risk that are easily updated using model parameter estimates from the existing literature. That is, developing models of risk that can be parameterised without raw data. This will require parameter estimates from different studies to be combined into one model, and therefore a common scale for risk is required. Two important scales in genetic association studies are the logit and the liability scale. We therefore have derived approximations to translate univariate effect size estimates for risk variables from the logit to the liability scale, and vice versa. This allows researchers to select a common scale, and transform any parameter estimates to this scale if required.

In particular, we focus on transforming effect size estimates from the logit to the liability scale. The liability scale offers a simple way to combine univariate effect size estimates for multiple risk factors into a joint effects model, with effect sizes not requiring re-calculation as more risk factors are identified and included in the model (assuming risk factors are independent). The same is not true for the logistic model. 

In the example R-script 'scale_transformation_code_HUMAN_HEREDITY_PAPER' we demonstrate how to transform effect size estimates from the logit scale to the liability scale. This is done using schizophrenia as an example. We find effect size estimates for 5 environmental risk factors on the liability scale using odds ratios reported in the literature. We create an environmental risk score (ERS) using these transformed effect size estimates and then estimate the risk of schizophrenia using this ERS and an existing schizophrenia PRS.

The R-script 'preprint_version_tss_paper_schizophrenia_eg' relates to the preprint version of the paper, and code does not match figure numbers in the final human heredity paper, and does not include an extra logistic model exploration.

