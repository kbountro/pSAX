# pSAX
Algorithm and scripts to implement the Kernel-based Probabilistic SAX method

# Table of Contents
1. [Introduction](#introduction)
2. [Modules Description](#files)
3. [Datasets](#datasets)
4. [Installation and Execution Instructions](#execution)


## Introduction <a name="introduction"></a>
The pSAX (Kernel-based Probabilistic SAX) [[1]](#1) method is an extension of the well-known SAX [[2]](#2) (Symbolic Aggregate Approximation) for time-series dimensionality reduction. The main contribution of the method is a SAX-based representation that adapts directly to the underlying probability distribution of the time-series data, thus providing a more accurate symbolic approximation. The accuracy has been measured and compared to the conventional SAX with the (significant for the indexing community) Tightness of Lower Bound metric, and also with the Mean Squared Error.


## Files Description <a name="files"></a>
This project consists of the following components:

* **demo:** Demo script (check [below](#execution) how to use). The code performs Monte-Carlo experiments for pSAX and SAX and then plots their Tightness of Lower Bound and Mean Squared Error. Each iteration consists of a comparison between a random "test" (target) subsequence from the dataset and another random "query" subsequence from the dataset.
* **tsPAA:** Time-series to PAA approximation. The PAA segments are single points (that is, the output in not really "segments", but the the values of the segments).
* **timeseries2symbol:** (c) 2003, Eamonn Keogh, Jessica Lin, Stefano Lonardi, Pranav Patel, Li Wei. Computes SAX representation of the data. The output are integer numbers, but should be seen as "symbols", not numbers.
* **min_paa_dist:**  Computes the lower-bounding distance, as defined in [[2]](#2).
* **plot_SAX:** Plots the SAX sequence against the raw time-series. Each symbol is plotted as a segment, with length equal to the PAA segment it was used. Also plots the estimated density function of the data and the "breakpoints" that are used to assign symbols to the PAA segments.
* **mvksdensity, statskcompute, statskernelinfo:** These are MATLAB's files (c) 2015-2016 The MathWorks, Inc. They are called from the built-in function 'ksdensity'. We tweaked them to i) allow to train from arbitrarily large number of samples (it was limited to 100 samples before) and ii) to fix the optimal smoothness parameter estimation for the Epanechnikov kernel, as it was set for the Gaussian kernel only. See https://www.mathworks.com/help/stats/ksdensity.html for more info.
* **lloydmax:** Lloyd-Max quantizer. Quantize according to a probability density function.
* **k-means++:** The k-means++ algorithm for initialization of k-means. Taken from the k-means file of Laurent S.: (https://www.mathworks.com/matlabcentral/fileexchange/28804-k-means), version 1.7.0.0


## Datasets <a name="datasets"></a>
We tested our method using some of the datasets available in https://www.cs.ucr.edu/~eamonn/iSAX/iSAX.html


## Installation and Execution Instructions <a name="execution"></a>
1. Download the project's source files.
2. Export as they are to a single folder.
3. Open MATLAB and load any 1D dataset. You should name the data variable as 'data'.
4. Call demo.m


## References
<a id="1">[1]</a> 
K. Bountrogiannis, G. Tzagkarakis and P. Tsakalides, "Data-driven Kernel-based Probabilistic SAX for Time Series Dimensionality Reduction," 2020 28th European Signal Processing Conference (EUSIPCO), Amsterdam, 2021, pp. 2343-2347, doi: 10.23919/Eusipco47968.2020.9287311.

<a id="2">[2]</a> 
J. Lin et al., “Experiencing SAX: A novel symbolic representation of time series”, Data Min. Knowl. Disc., vol. 15, no. 2, pp. 107–144, 2007

## LICENSE
This source code can be used for non-commercial purposes only. Its utilization must acknowledge and cite the following publication:
K. Bountrogiannis, G. Tzagkarakis and P. Tsakalides, "Data-driven Kernel-based Probabilistic SAX for Time Series Dimensionality Reduction," 2020 28th European Signal Processing Conference (EUSIPCO), Amsterdam, 2021, pp. 2343-2347, doi: 10.23919/Eusipco47968.2020.9287311. 
