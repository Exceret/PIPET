# PIPET

Phenotypic information based on bulk data predicts relevant subpopulations in single cell data

Please use [SigBridgeR](https://github.com/WangLabCSU/SigBridgeR) for the fork version of PIPET. Bug reports and feature requests are welcomed at [SigBridgeR-issues](https://github.com/WangLabCSU/SigBridgeR/issues).


## Introduction

`PIPET` can be used to predict relevant subpopulations in single-cell data from phenotypic information in bulk data. You can use known feature vectors of phenotypic information to preform PIPET() or create feature vectors via the function in the `PIPET`. This package also provides commonly used downstream analysis for creating customizable visualization results. The workflow of `PIPET` is shown in the following Figure:

<p align="center">
<img src=Figure_PIPET.jpg height="900" width="640">
</p>

## Installation

To install this package, start R (*version >= 4.3.1*) and enter:

```R
if (!requireNamespace("pak")) {
  install.packages("pak")
}
pak::pkg_install("Exceret/PIPET")
```



