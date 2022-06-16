% pSAX representation, as introduced and studied in the following papers:
%
%   [1] K. Bountrogiannis, G. Tzagkarakis and P. Tsakalides, 
%   "Data-driven Kernel-based Probabilistic SAX for Time Series Dimensionality Reduction,"
%   2020 28th European Signal Processing Conference (EUSIPCO), 2021.
%
%   [2] K. Bountrogiannis, G. Tzagkarakis and P. Tsakalides, 
%   "Distribution  Agnostic Symbolic Representations for Time Series Dimensionality Reduction and Online Anomaly Detection," 
%   in IEEE Transactions on Knowledge and Data Engineering, 2022.
%
%   NOTE: This version scans the dataset and extracts all subsequences,
%   even overlapping ones. This version is basically used for discord detection in [2].
%   For non-overlapping windows, check the other file (pSAX.m).
%
%   The fundamental difference between pSAX.m and pSAX_overlap.m is
%   that pSAX.m treats the dataset as a whole, whereas pSAX_overlap.m
%   treats the dataset as a bunch of separate subsequences.
%
%   Author: Konstantinos Bountrogiannis
%   Contact: kbountrogiannis@gmail.com
%   Date: June 2022

function [str] = pSAX_overlap(data, training_len, win_size, paa_size, alphabet_size, normalize)

% Inputs:   data: data :)
%           training_len:   the number of samples used for training
%           win_size:   length of sliding window
%           paa_size:   the size (dimensionality) of the PAA approximation of
%                       the subsequence in the sliding window
%           alphabet_size:  quantizer's codebook size
%           normalize:  (0 or 1). Option to normalize each subsequence prior to
%                       processing. Normally, this is true (1).
%           
% Output:   str:    pSAX symbolic sequence


    % We will first create the training set. Note that we treat our dataset
    % as a bunch of seperate subsequences. Therefore, the normalization is
    % performed independently on each subsequence.
    training_PAA = zeros(training_len-win_size+1,paa_size);
    for i = 1 : training_len-win_size+1
        sub_section = data(i:i + win_size -1);
        if normalize
            if std(sub_section)<0.001
                sub_section = sub_section - mean(sub_section);
            else
                sub_section = (sub_section - mean(sub_section))/std(sub_section);
            end
        end
        training_PAA(i,:) = tsPAA(sub_section,paa_size);
    end
    
    % The training will be performed on the PAA sequence. We can
    % concatenate it.
    training_PAA = training_PAA(:)';

    % Kernel density estimation from the training set
    [f,x] = ksdensity(training_PAA,'npoints',training_len,'Kernel','epanechnikov');
    
    % Lloyd-Max quantization
    [~,init_codewords] = kmeanspp(training_PAA,alphabet_size); % Initialize codebook with k-means++
    init_codewords = sort(init_codewords,2);
    [~,cutlines] = lloydmax(f,x,alphabet_size,[],init_codewords);

    % Discretization as usual, with the custom cutlines
    str = timeseries2symbol(data, win_size, paa_size, cutlines, normalize, 1) - 1;
            
end