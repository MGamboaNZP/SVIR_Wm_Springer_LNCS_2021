# Moments of Wm in a Stochastic SVIR Model

This repository contains MATLAB implementations for the computation of the first and second moments of the random variable Wm in a stochastic SVIR epidemic model with imperfect vaccine and an external source of infection.

The code accompanies and supports the results presented in the following reference:

**Gamboa, M., López-García, M., & Lopez-Herrero, M. J. (2021).**  
*A stochastic SVIR model with imperfect vaccine and external source of infection.*  
In *European Workshop on Performance Engineering*, pp. 197–209.  
Springer International Publishing.  
https://doi.org/10.1007/978-3-030-91825-5_12

If you use this code in your research, please cite the following publication:

Gamboa, M., López-García, M., & Lopez-Herrero, M. J. (2021).
A stochastic SVIR model with imperfect vaccine and external source of infection.
European Workshop on Performance Engineering, 197–209.
https://doi.org/10.1007/978-3-030-91825-5_12

---

## Repository contents

The repository includes the following MATLAB code:

- `momentoW_M.m`  
  Computes the first and second moments of the warning vaccination level using recursive equations derived from the stochastic SVIR model.

The implementation closely follows the theoretical development in the cited paper and is intended for research and reproducibility purposes.

---

## Model description

The model considers a stochastic SVIR epidemic process with:
- Imperfect vaccination,
- External source of infection,


The algorithm computes recursively:
- The first moment (expected value),
- The second moment,

of the warning vaccination level, conditional on the initial state of the system.

---

## Requirements

- MATLAB (R2018a or later recommended)
- No additional toolboxes are required

---

## Usage example

The main function can be called as follows:

```matlab
[K1, K2] = momentoW_M(v0, s0, i0, m, N, gamma, xi, beta, h);
