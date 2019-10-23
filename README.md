[![DOI](https://zenodo.org/badge/89879581.svg)](https://zenodo.org/badge/latestdoi/89879581)

# w4mclassfilter

This is a Galaxy-oriented R package that can be used to filter Workflow4Metabolomics ("W4M", Giacomoni *et al.*, 2014) data matrix and metadata by sample-class.
(The data matrix and metadata files are written by the W4M wrapper of XCMS, Smith *et al.*, 2006.)
For a brief introduction to this package, please see the vignette: (https://github.com/HegemanLab/w4mclassfilter/blob/master/vignettes/w4mclassfilter.Rmd).

This package has been "wrapped" as a Galaxy tool
  - the "wrapping project" is here: (https://github.com/HegemanLab/w4mclassfilter_galaxy_wrapper)
  - the "w4mclassfilter" tool is in the toolshed at (https://toolshed.g2.bx.psu.edu/repository?repository_id=5f24951d82ab40fa)

## Installing the w4mclassfilter R package

To use this package in R, download the tar.gz file for the desired version from (https://github.com/HegemanLab/w4mclassfilter/releases), and install it into R one of the following ways:

### From within RStudio, 

```
From the "Tools" menu, choose "Install packages..."
Change "Install from:" to "Package Archive File (.tar.gz)", and a "file picker" dialog will pop up.
Navigate to and select w4mclassfilter_blahblah.tar.gz
Click the 'Install' button.
```

### From within R,

```r
install.packages("path/to/w4mclassfilter_blahblah.tar.gz", repos = NULL, type="source")
```

### From the command line,

```bash
cd directory_that_contains_w4mclassfilter_blahblah.tar.gz
R CMD INSTALL w4mclassfilter_blahblah.tar.gz
```

## License

MIT License

Copyright (c) 2017 Arthur C Eschenlauer

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## NEWS

### CHANGES IN VERSION 0.98.14

* Built with R 3.6.1 (2019-07-05) "Action of the Toes"

#### NEW FEATURES

* Enhancement https://github.com/HegemanLab/w4mclassfilter_galaxy_wrapper/issues/6 - "Provide sort options for features and samples"

#### INTERNAL MODIFICATIONS

* Support the enhancement

### CHANGES IN VERSION 0.98.13

* Built with R 3.6.1 (2019-07-05) "Action of the Toes"

#### NEW FEATURES

* none

#### INTERNAL MODIFICATIONS

* Fix https://github.com/HegemanLab/w4mclassfilter/issues/5 - "MAXFEAT does not work for imputation functions that include data transformation"


### CHANGES IN VERSION 0.98.12

#### NEW FEATURES

* Enhancement https://github.com/HegemanLab/w4mclassfilter/issues/4 - "add and test no-imputation and centering-imputation functions":
  - Support no imputation.
  - Support imputating missing feature-intensities as median intensity for the corresponding feature.  

#### INTERNAL MODIFICATIONS

* Support and tests for new features.
* Built with R 3.6.1 (2019-07-05) "Action of the Toes"
* Version numbers 0.98.10 and 0.98.10 were skipped to synchronize with w4mclassfilter_galaxy_wrapper


### CHANGES IN VERSION 0.98.9

#### NEW FEATURES

* none

#### INTERNAL MODIFICATIONS

* Fix https://github.com/HegemanLab/w4mclassfilter/issues/3 - "a single column is insufficient to calculate a variance"
* Built with R version 3.4.1 (2017-06-30) "Single Candle"

### CHANGES IN VERSION 0.98.8

#### NEW FEATURES

* none

#### INTERNAL MODIFICATIONS

* Fix https://github.com/HegemanLab/w4mclassfilter/issues/1 - "bad error handling when input file cannot be found"
* Add tests for missing and invalid input files
* Built with R version 3.4.1 (2017-06-30) "Single Candle"

### CHANGES IN VERSION 0.98.7

#### NEW FEATURES

* First column of sample metadata is by default renamed to "sampleMetadata" unless 
  argument 'name_smplmetadata_col1' is supplied and set to FALSE.

#### INTERNAL MODIFICATIONS

* none

### CHANGES IN VERSION 0.98.6

#### NEW FEATURES

* Add field filters for variableMetadata as csv; default empty.
  - E.g. for mz &gt; 200 and mz &lt; 800 and rt &gt; 150 and rt &lt; 1100,
    `mz:200:800,rt:150:1100`

#### INTERNAL MODIFICATIONS

* Sort dataMatrix to match order of sampleMetadata and variableMetadata


### CHANGES IN VERSION 0.98.3

#### NEW FEATURES

* Galaxy-support unchanged.
* Provide for more flexible input and output from env, list, or matrix/data.frame. 

#### INTERNAL MODIFICATIONS

* Support and tests for new features.


### CHANGES IN VERSION 0.98.2

#### NEW FEATURES

* Added support for R-flavored regular expression pattern matching.
* Empty classes argument or zero-length class\_column result in no samples filtered out. 

#### INTERNAL MODIFICATIONS

* Support and tests for new features.

#### BUILD NOTES

* This package was built from ~/src/w4mclassfilter with R 3.3.1 under bash and miniconda as follows:

```
  $ ~/miniconda2/bin/conda create -n r3.3.1 r-base=3.3.1
  $ source ~/miniconda2/bin/activate r3.3.1
  (r3.3.1) $ R
  > library(devtools)
  > document()
  > test()
  > check()
  > build()
  > install()

```

Result is at ~/src/w4mclassfilter\_0.98.2.tar.gz

### CHANGES IN VERSION 0.98.1

#### NEW FEATURES

* First release - Wrap the w4mclassfilter R package that implements filtering of W4M data matrix, variable metadata, and sample metadata by class of sample.
* *dataMatrix* *is* modified by the tool, so it *does* appear as an output file
* *sampleMetadata* *is* modified by the tool, so it *does* appear as an output file
* *variableMetadata* *is* modified by the tool, so it *does* appear as an output file

#### INTERNAL MODIFICATIONS

* none

## Citations

Smith, Colin A. and Want, Elizabeth J. and O’Maille, Grace and Abagyan, Ruben and Siuzdak, Gary (2006). XCMS: Processing Mass Spectrometry Data for Metabolite Profiling Using Nonlinear Peak Alignment, Matching, and Identification. In Analytical Chemistry, 78 (3), pp. 779–787. [doi:10.1021/ac051437y](http://dx.doi.org/10.1021/ac051437y)

Giacomoni, F. and Le Corguille, G. and Monsoor, M. and Landi, M. and Pericard, P. and Petera, M. and Duperier, C. and Tremblay-Franco, M. and Martin, J.-F. and Jacob, D. and et al. (2014). Workflow4Metabolomics: a collaborative research infrastructure for computational metabolomics. In Bioinformatics, 31 (9), pp. 1493–1495. [doi:10.1093/bioinformatics/btu813](http://dx.doi.org/10.1093/bioinformatics/btu813)
