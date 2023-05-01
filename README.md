# Panel Local Projection

This repository hosts code that accompanies

* Ziwei Mei, Liugang Sheng, Zhentao Shi (2023), "[Implicit Nickell Bias in Panel Local Projection](https://arxiv.org/abs/2302.13455)", _arxiv:2302.13455_. 

It offers an `R` function to implement the panel local projection that includes two methods: 

- The conventional _fixed effect_ (FE) estimator. This is the default estimator in most applied works, but we find that it incurs Nickell bias
- The _half-panel jackknife_ (HJ) estimator eliminates the asymptotical bias and delivers valid statistical inference.



### Contents  

- `LP_panel.R` includes the main function `LP_panel()` with detailed instructions. 
- `example.R` provides a running example to demonstrate the usage of the function `LP_panel()`.



### Contributors 

Ziwei Mei, Zhentao Shi, and Shen Shu
