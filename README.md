# Code for manuscript: "Scarce Data, Noisy Inferences, and Overfitting: The Hidden Flaws in Ecological Dynamics Modeling" 

## Python codes
Codes that implement model reduction for the generalized lotka volterra model

### Installation
- Create and environment
```python
cd python
python -m venv mbam_venv
```
Activate it
```python
source mbam_venv/bin/activate
```
Upgrade pip, just in case
```python
pip install --upgrade pip
```
Install required packages
```python
pip install -r requirements.txt
```
Run one example
```python
python fourpop_mbam_reduction.py
```

## R codes

We simulate the deterministic generalized Lotka-Volterra model using the R library `deSolve`. For inference, the state-of-art-bayesian engine, `stan`.

All the auxiliary functions are packed in file `R/lotka_volterra_stan_functions.R`. The file `R/batch_434.R` illustrates how to choose a seed (434 in this case), noise levels and population sizes to reproduce the figures in the article. This script calls `R/lotka_volterra_rk4_stan.R` that makes the bayesian inference.

Required libraries: `deSolve`, `rstan`, `ggplot2`, `psych`, `deSolve`, `rstan`, `ggplot2`, `psych`, `bayesplot`. 


