# Kernel-based probabilistic SAX (pSAX)
MATLAB implementation of the pSAX time-series symbolic representation. Tested with MATLAB R2019a.

# Table of Contents
1. [Introduction](#introduction)
2. [Modules Description](#files)
3. [Datasets](#datasets)
4. [Installation and Execution Instructions](#execution)


## Introduction <a name="introduction"></a>
The pSAX (Kernel-based Probabilistic SAX) [[1]](#1), [[2]](#2) method is an extension of the well-known SAX [[2]](#3) (Symbolic Aggregate Approximation) for time-series dimensionality reduction. The main contribution of the method is a SAX-based representation that adapts directly to the underlying probability distribution of the time-series data, thus providing a more accurate symbolic approximation. The accuracy has been measured and compared to the conventional SAX with the (significant for databases performance) Tightness of Lower Bound metric, and also with the Mean Squared Error.


## Files Description <a name="files"></a>
This project consists of the following components:

* **pSAX:** Main function.
* **tsPAA:** (c) 2003, Eamonn Keogh, Jessica Lin, Stefano Lonardi, Pranav Patel, Li Wei. Time-series to PAA approximation.
* **timeseries2symbol:** (c) 2003, Eamonn Keogh, Jessica Lin, Stefano Lonardi, Pranav Patel, Li Wei. Computes SAX representation of the data. The output are integer numbers.
* **map_to_string:** Discretization of PAA sequence with the custom quantization intervals.
* **mvksdensity, statskcompute, statskernelinfo:** These are MATLAB's source files, (c) 2015-2016 The MathWorks, Inc. They are called from the built-in function 'ksdensity'. We tweaked them to i) allow to estimate arbitrarily large number of density points (it was limited to 100 before) and ii) to fix the optimal smoothness parameter estimation for the Epanechnikov kernel, as it was set for the Gaussian kernel only. See https://www.mathworks.com/help/stats/ksdensity.html for more info.
* **lloydmax:** Lloyd-Max quantizer. Quantize according to a probability density function.
* **k-means++:** The k-means++ algorithm for initialization of k-means. Taken from the k-means file of Laurent S.: (https://www.mathworks.com/matlabcentral/fileexchange/28804-k-means), version 1.7.0.0


## Datasets <a name="datasets"></a>
A large collection of datasets is available in https://www.cs.ucr.edu/~eamonn/iSAX/iSAX.html


## Installation and Execution Instructions <a name="execution"></a>
1. Download the project's source files.
2. Export as they are to a single folder.
3. Call pSAX with the appropriate inputs.


## References
<a id="1">[1]</a> 
K. Bountrogiannis, G. Tzagkarakis and P. Tsakalides, "Data-driven Kernel-based Probabilistic SAX for Time Series Dimensionality Reduction," 2020 28th European Signal Processing Conference (EUSIPCO), Amsterdam, pp. 2343-2347, 2021, doi: 10.23919/Eusipco47968.2020.9287311.

<a id="2">[2]</a> 
K. Bountrogiannis, G. Tzagkarakis and P. Tsakalides, "Distribution Agnostic Symbolic Representations for Time Series Dimensionality Reduction and Online Anomaly Detection," in IEEE Transactions on Knowledge and Data Engineering, doi: 10.1109/TKDE.2022.3174630.

<a id="3">[3]</a> 
J. Lin et al., “Experiencing SAX: A novel symbolic representation of time series”, Data Min. Knowl. Disc., vol. 15, no. 2, pp. 107–144, 2007, doi: 10.1007/s10618-007-0064-z

## License
**This code is released under GPL v.3.0. If you use this code for academic works, please cite at least one of the following publications:**

K. Bountrogiannis, G. Tzagkarakis and P. Tsakalides, "Data-driven Kernel-based Probabilistic SAX for Time Series Dimensionality Reduction," 2020 28th European Signal Processing Conference (EUSIPCO), Amsterdam, 2021, pp. 2343-2347, doi: 10.23919/Eusipco47968.2020.9287311.

K. Bountrogiannis, G. Tzagkarakis and P. Tsakalides, "Distribution Agnostic Symbolic Representations for Time Series Dimensionality Reduction and Online Anomaly Detection," in IEEE Transactions on Knowledge and Data Engineering, doi: 10.1109/TKDE.2022.3174630.
