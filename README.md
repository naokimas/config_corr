# config_corr
Configuration model for correlation/covariance matrices

This package contains two algorithms.
* Altogihm DMCC solves a convex optimization problem and usually better than the other naive gradient descent algorithm. Therefore, we recommend DMCC over the naive gradient descent algorithm.
* DMCC, including the performance test, is documented in the accompanying pdf file.
* DMCC is provided in Python only.
* The naive gradient descent algorithm is provided in Python and MATLAB.

## DMCC algorithm in Python

```python
python3 test_dmcc.py
```

## Naive gradient descent algorithm in Python

```python
python3 test_naive_gradient_descent.py
```


## Naive gradient descent algorithm in MATLAB

```MATLAB
test_naive_gradient_descent
```

* test_naive_gradient_descent.m calls max_ent_config_naive_gradient_descent.m
* The data set to be used should be specified within test_naive_gradient_descent.m
* Variable curr_dir should be modified according to where you place the data set.
