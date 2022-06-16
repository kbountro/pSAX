% pSAX representation, as introduced and studied in the following papers:
%
%   K. Bountrogiannis, G. Tzagkarakis and P. Tsakalides, 
%   "Data-driven Kernel-based Probabilistic SAX for Time Series Dimensionality Reduction,"
%   2020 28th European Signal Processing Conference (EUSIPCO), 2021.
%
%   K. Bountrogiannis, G. Tzagkarakis and P. Tsakalides, 
%   "Distribution  Agnostic Symbolic Representations for Time Series Dimensionality Reduction and Online Anomaly Detection," 
%   in IEEE Transactions on Knowledge and Data Engineering, 2022.

function [str] = pSAX_overlap(data, training_len, win_size, paa_size, alphabet_size, normalize)

% Inputs:   data: data :)
%           training_len:   the number of samples used for density estimation
%           win_size:   length of sliding window
%           paa_size:   the size (dimensionality) of the PAA approximation of
%                       the subsequence in the sliding window
%           alphabet_size:  quantizer's codebook size
%           normalize:  (0 or 1). Option to normalize each subsequence prior to
%                       processing. Normally, this is true (1).
%           
% Output:   str:    pSAX symbolic sequence
%
% Author: Konstantinos Bountrogiannis
% Contact: kbountrogiannis@gmail.com
% Date: June 2022

    % For pSAX we need a training set. Here, we opt to train on the first
    % section of our data.
    paa_train = zeros(training_len-win_size+1,paa_size);
    for i = 1 : training_len-win_size+1
        sub_section = data(i:i + win_size -1);
        if normalize
            sub_section_mean = mean(sub_section);
            sub_section_std = std(sub_section);
            if sub_section_std>0.001
                sub_section = (sub_section - sub_section_mean)/sub_section_std;
            else
                sub_section = sub_section - sub_section_mean;
            end
        end
        paa_train(i,:) = tsPAA(sub_section,paa_size);
    end
    paa_train = paa_train(:)';

    % Kernel density estimation from the training array
    [f,x] = ksdensity(paa_train,'npoints',training_len,'Kernel','epanechnikov');
    
    % Lloyd-Max quantization
    [~,init_codewords] = kmeanspp(paa_train,alphabet_size); % Initialize codebook with k-means++
    init_codewords = sort(init_codewords,2);
    [~,cutlines] = lloydmax(f,x,alphabet_size,[],init_codewords);
    % now 'cutlines' are the trained boundaries of the discretization
    % intervals

    % Discretization as usual, with the custom cutlines
    str = timeseries2symbol(data, win_size, paa_size, cutlines, normalize, 1) - 1;
            
end

