# Moments of Wm in a Stochastic SVIR Model

This repository contains R scripts associated with the publication:

Gamboa, M., López-García, M., & Lopez-Herrero, M. J. (2021).
A stochastic SVIR model with imperfect vaccine and external source of infection.
European Workshop on Performance Engineering, 197–209.
https://doi.org/10.1007/978-3-030-91825-5_12

---

## Repository contents

The repository includes the following MATLAB code:

- `momentoW_M.m`  
  Computes moments of the random variable Wm in a stochastic SVIR model.
The implementation closely follows the theoretical development in the cited paper and is intended for research and reproducibility purposes.

---



## Usage example

The main function can be called as follows:

```matlab
[K1, K2] = momentoW_M(v0, s0, i0, m, N, gamma, xi, beta, h);
