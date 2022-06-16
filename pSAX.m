%   pSAX representation, as introduced and studied in the following papers:
%
%   [1] K. Bountrogiannis, G. Tzagkarakis and P. Tsakalides, 
%   "Data-driven Kernel-based Probabilistic SAX for Time Series Dimensionality Reduction,"
%   2020 28th European Signal Processing Conference (EUSIPCO), 2021.
%
%   [2] K. Bountrogiannis, G. Tzagkarakis and P. Tsakalides, 
%   "Distribution  Agnostic Symbolic Representations for Time Series Dimensionality Reduction and Online Anomaly Detection," 
%   in IEEE Transactions on Knowledge and Data Engineering, 2022.
%
%   NOTE: the dataset is scanned with non-overlapping windows. If you want to 
%   scan the dataset and extract all subsequences, even overlapping ones,
%   check the other file (pSAX_overlapping). 
%
%   The fundamental difference between pSAX.m and pSAX_overlap.m is
%   that pSAX.m treats the dataset as a whole, whereas pSAX_overlap.m
%   treats the dataset as a bunch of separate subsequences.
%
%   Author: Konstantinos Bountrogiannis
%   Contact: kbountrogiannis@gmail.com
%   Date: June 2022

function [str] = pSAX(data, training_len, dim_ratio, alphabet_size, normalize)

% Inputs:   data: data :)
%           training_len:   the number of samples used for density estimation
%           dim_ratio:  dimensionality reduction ratio. For example, if
%                       dim_ratio is 1/3, then every 3 samples of the
%                       dataset will be mapped to 1 symbol in the symbolic
%                       sequence.
%           alphabet_size:  quantizer's codebook size
%           normalize:  (0 or 1). Option to normalize each subsequence prior to
%                       processing. Normally, this is true (1).
%           
% Output:   str:    pSAX symbolic sequence.
%

    % First, remove the remainder of the dataset, such that it fits exactly.
    data_len = length(data);
    data_nseg = floor(dim_ratio*data_len);
    data_len = data_len - mod(data_len,data_nseg);
    data = data(1:data_len);
    % Do the same for the training set.
    training_nseg = floor(dim_ratio*training_len);
    training_len = training_len - mod(training_len,training_nseg);
    training_set = data(1:training_len); % Note that we opt to train on the first section of our data.

    % We will first normalize our data. Notice that we treat our dataset
    % (and the training set) as a whole, not as seperate subsequences.
    if normalize
        if std(training_set) <0.001
            data = data - mean(data);
            training_set = training_set - mean(training-set);
        else
            data = (data-mean(data))/std(data);
            training_set = (training_set-mean(training_set))/std(training_set);
        end
    end
    
    % The training will be performed on the PAA sequence
    training_PAA = tsPAA(training_set,training_nseg);

    % Kernel density estimation from the training set
    [f,x] = ksdensity(training_PAA,'npoints',training_len,'Kernel','epanechnikov');
    
    % Lloyd-Max quantization
    [~,init_codewords] = kmeanspp(training_PAA,alphabet_size); % Initialize codebook with k-means++
    init_codewords = sort(init_codewords,2);
    [~,cutlines] = lloydmax(f,x,alphabet_size,[],init_codewords);
    
    % Discretization as usual, but with the learned cutlines
    str = timeseries2symbol(data, data_len, data_nseg, cutlines, normalize, 1) - 1;
    str = reshape(str,[1,data_nseg]);
            
end

