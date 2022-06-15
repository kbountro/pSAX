% Copyright and terms of use (DO NOT REMOVE):
% The code is made freely available for non-commercial uses only, provided that the copyright 
% header in each file not be removed, and suitable citation(s) (see below) be made for papers 
% published based on the code.
%
% The code is not optimized for speed, and we are not responsible for any errors that might
% occur in the code.
%
% The copyright of the code is retained by the authors.  By downloading/using this code you
% agree to all the terms stated above.
%
%   Lin, J., Keogh, E., Lonardi, S. & Chiu, B. 
%   "A Symbolic Representation of Time Series, with Implications for Streaming Algorithms." 
%   In proceedings of the 8th ACM SIGMOD Workshop on Research Issues in Data Mining and 
%   Knowledge Discovery. San Diego, CA. June 13, 2003. 
%
%
%   Lin, J., Keogh, E., Patel, P. & Lonardi, S. 
%   "Finding Motifs in Time Series". In proceedings of the 2nd Workshop on Temporal Data Mining, 
%   at the 8th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining. 
%   Edmonton, Alberta, Canada. July 23-26, 2002
%
% Copyright (c) 2003, Eamonn Keogh, Jessica Lin, Stefano Lonardi, Pranav Patel, Li Wei. All rights reserved.
%

function PAA = tsPAA(data,n)
% Compute the Piecewise Aggregate Approximation (PAA) of a Time-Series
% The output is stretched to match the size of the input
%
% Input:    data : Time-Series input
%           n    : number of segments to split the input
% Output:   PAA  : Computed PAA

N = length(data);
win_size = floor(N/n);
d = min(size(data));

% Take care of the special case where there is no dimensionality reduction
if N == n
    PAA = data;

else
        % N is not dividable by n
        if (N/n - floor(N/n)) %#ok<*BDLOG>
            temp = zeros(n, N);
            for j = 1 : n
                temp(j, :) = data;
            end
            expanded_data = reshape(temp, 1, N*n);
            PAA = mean(reshape(expanded_data, N, n));
        % N is dividable by n
        else
            PAA_temp = reshape(data,win_size,n);
            PAA = mean(PAA_temp);
    %         [~,maxidx] = max(abs(PAA),[],'linear');
    %         PAA = PAA(maxidx);
        end
end

end

