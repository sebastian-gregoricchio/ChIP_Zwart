---
title: "changeLog"
---

# ChIP_Zwart
[![forks](https://img.shields.io/github/forks/sebastian-gregoricchio/ChIP_Zwart?style=social)](https://github.com/sebastian-gregoricchio/ChIP_Zwart/fork)


#### [v0.1.3](https://github.com/sebastian-gregoricchio/ChIP_Zwart/releases/tag/0.1.3) - October 24<sup>th</sup> 2024
* Bug fixing for pattern detection of sample names in the peak calling analyses
* Added feature to get the counts from multiBigWigSummary also in tab format


#### [v0.1.2](https://github.com/sebastian-gregoricchio/ChIP_Zwart/releases/tag/0.1.2) - September 7<sup>th</sup> 2023
* Optimization of the threads management in the workflow snakefiles
* Adding a comment for eventual bug fixing for the path to use to load python modules in the snakemake (workflow) files
* Added `$CONDA_PREFIX/bin/` behind each tool to avoid wrong assignment of the tool path


#### [v0.1.1](https://github.com/sebastian-gregoricchio/ChIP_Zwart/releases/tag/0.1.1) - March 31<sup>st</sup> 2023
* Updated conda environment yaml file to avoid to point to the wrong python/R libs: `conda-forge` package `conda-ecosystem-user-package-isolation=1.0=ha770c72_1` added


#### [v0.1.0](https://github.com/sebastian-gregoricchio/ChIP_Zwart/releases/tag/0.1.0) - March 29<sup>th</sup> 2023
* First release
