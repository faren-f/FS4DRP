# Comparative evaluation of knowledge-based feature reduction methods for drug response prediction

Personalized medicine aims to tailor medical treatments to individual patients, and predicting drug responses from molecular profiles using machine learning (ML) is crucial for this goal. However, the high dimensionality of molecular profiles compared to the limited number of samples presents significant challenges. Knowledge-based feature selection methods are particularly suitable for drug response prediction as they leverage biological insights to reduce dimensionality and enhance model interpretability. This study presents the first comparative evaluation of commonly used knowledge-based feature reduction methods using cell line and tumor data. Our findings indicate that transcription factor activities outperform other methods in predicting drug responses, effectively distinguishing between sensitive and resistant tumors for seven out of 20 drugs evaluated.


## Instructions

### Install packages:

In *R*
```
glmnet(4.1.8)
caret(6.0.94)
keras(2.15.0)
tensorflow(2.16.0)
randomForest(4.7.1.1)
progeny(1.24.0)
reactome.db(1.86.2)
sva(3.50.0)
caTools(1.18.2)
rtracklayer(1.62.0)
readxl(1.4.3)
parallel(4.3.2)
```
In *Python*

```
pandas
decoupler
```
*hint:* While the main pipeline is implemented in *R*, the *Python* packages are required to implement transcription factor activists only.

### Perform analysis:

In this work, we compare four different feature reduction methods including

1. Landmark genes
2. Drug pathway genes
3. Pathway activities
4. TF activities

together with five ML models including

1. Ridge regression
2. Lasso regression
3. Elastic net
4. Random forest
5. Multilayer perceptron

Our analysis consists of two scenarios:

**1. Cross-validation on cell lines.** In this scenario, train and test sets for ML models are obtained from cell line data (PRISM dataset) using cross-validation. The script names for this analysis follow the pattern ```PRISM_{feature_selection_method}.R```.

**2. Validation on tumour samples.** In this scenario, train and test sets for ML models are obtained from cell line data (PRISM dataset) and tumor data (TCGA dataset), respectively. The script names for this analysis follow the pattern ```TrainPRISM_TestTCGA_{feature_selection_method}.R```.


*hint:* The above mentioned scripts can be found in the ```main``` folder.

In the following, the step-by-step instructions to run the pipline and obtain the results are described:

1. *Download raw data.* see [here](data/raw_data/README.md)
2. *data_preprocessing:*




