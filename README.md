# pSAX
pSAX necessary files and demo script

# Table of Contents
1. [Introduction](#introduction)
2. [Modules Description](#files)
3. [Datasets](#datasets)
4. [Installation Instructions](#execution)

## Introduction <a name="introduction"></a>
The pSAX (Kernel-based Probabilistic SAX) method is an extension of the well-known SAX (Symbolic Aggregate Approximation) for time-series dimensionality reduction. The main contribution of the method is a SAX-based representation that adapts directly to the underlying probability distribution of the time-series data, thus providing a more accurate symbolic approximation. The accuracy has been measured and compared to the conventional SAX with the (significant for the indexing community) Tightness of Lower Bound metric, but also with the Mean Squared Error.


## Files Description <a name="files"></a>
This project consists of the following components:

* **Tensor:** Class for the construction and management of 3D data structures (the tensors).
* **TMacPar:** Class for reconstructing low-rank tensor measurements using the TMac algorithm.
* **Matlab:** Class that provides useful mathematical functions, operations in three-dimensional structures, as well as type conversion between its various mathematical libraries.
* **IO:**  Class for reading and writing .csv files. Methods of this class are used to record the reconstruction results in the res\ folder, for post-processing evaluation purposes.
* **Main:** Class that includes the functions for executing TC in predefined data, as well as a function for user-defined data reconstruction.


## Datasets <a name="datasets"></a>
We tested our method using some of the datasets available in https://www.cs.ucr.edu/~eamonn/iSAX/iSAX.html


## Installation Instructions <a name="execution"></a>
1. Download the project's source files.
2. Export as they are to a single folder.
3. Call sax_demo.m with MATLAB (it has been tested with versions 2018b and later).

