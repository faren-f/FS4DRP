# A comparative analysis of feature reduction methods for drug response prediction

Personalized medicine aims to tailor medical treatments to individual patients, and predicting drug responses from molecular profiles using machine learning (ML) is crucial for this goal. However, the high dimensionality of molecular profiles compared to the limited number of samples presents significant challenges. Knowledge-based feature selection methods are particularly suitable for drug response prediction as they leverage biological insights to reduce dimensionality and enhance model interpretability. This study presents the first comparative evaluation of commonly used knowledge-based feature reduction methods using cell line and tumor data. Our findings indicate that transcription factor activities outperform other methods in predicting drug responses, effectively distinguishing between sensitive and resistant tumors for seven out of 20 drugs evaluated.


## Instructions

### Install packages:

In *R*
```


```
In *Python*

```

```

### Run pipeline:

In this work, we compare four different feature selection methods including
1. 
2.
3.
4.

together with five ML models including
1.
2.
3.
4.
5.

Our analysis consists of two scenarios:

**1. Cross-validation on cell lines.** In this scenario, train and test sets for ML models are obtained from cell line data (PRISM dataset) using cross-validation. The script names for this analysis follow the pattern ```PRISM_{feature_selection_method}```.

**2. Validation on tumour samples.**

