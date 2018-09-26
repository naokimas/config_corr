# config_corr
Configuration model for correlation/covariance matrices

This package contains two algorithms.
* Altogihm DMCC solves a convex optimization problem and usually better than the other naive gradient descent algorithm. Therefore, we recommend DMCC over the naive gradient descent algorithm.
* DMCC, including the performance test, is documented in the accompanying [DMCC-algorithm.pdf](https://github.com/naokimas/config_corr/blob/master/DMCC-algorithm.pdf).
* DMCC is provided in Python only.
* The naive gradient descent algorithm is provided in Python and MATLAB.

## Install from source

configcorr has the following dependencies:

Python >= 3.4  
setuptools >= 1.4  
NumPy >= 1.8  
CVMPY

configcorr automatically installs NumPy and CVMPY.

Clone the configcorr git repository by running

```python
git clone https://github.com/naokimas/config_corr.git
```

Then, navigate to the top-level of the cloned directory and run

```python
python setup.py install
```

## DMCC algorithm in Python

```python
python3 test_dmcc.py
```

* test_dmcc.py calls max_ent_config_dmcc.py

## Naive gradient descent algorithm in Python

```python
python3 test_naive_gradient_descent.py
```

* test_naive_gradient_descent.py calls max_ent_config_naive_gradient_descent.py


## Naive gradient descent algorithm in MATLAB

```MATLAB
test_naive_gradient_descent
```

* test_naive_gradient_descent.m calls max_ent_config_naive_gradient_descent.m
* The data set to be used should be specified within test_naive_gradient_descent.m
* Variable curr_dir should be modified according to where you place the data set.
