# Panel Local Projection


This repository hosts an R package that accompanies

* Ziwei Mei, Liugang Sheng, Zhentao Shi (2023), "[Nickell Bias in Panel Local Projection](https://arxiv.org/abs/2302.13455)", _arxiv:2302.13455_. 

It offers an `R` function to implement the panel local projection that includes two methods: 

- The _fixed effect_ (FE) estimator. This is the default estimator in most applied works, but we find that it incurs Nickell bias
- The _split-panel jackknife_ (SPJ) estimator eliminates the asymptotical bias and delivers valid statistical inference.

```
devtools::install_github("zhentaoshi/panel-local-projection")
library("pLP")
```

### Contributors 

Ziwei Mei, Zhentao Shi, and Shen Shu


### License

This work is licensed under
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
