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
% This function takes in a time series and convert it to string(s).
% There are two options:
%   1. Convert the entire time series to ONE string
%   2. Use sliding windows, extract the subsequences and convert these subsequences to strings
%
% For the first option, simply enter the length of the time series as "N"
%   ex. We have a time series of length 32 and we want to convert it to a 8-symbol string,
%       with alphabet size 3:
%       timeseries2symbol(data, 32, 8, 3)
% For the second option, enter the desired sliding window length as "N"
%   ex. We have a time series of length 32 and we want to extract subsequences of length 16 using
%       sliding windows, and convert the subsequences to 8-symbol strings, with alphabet size 3:
%       timeseries2symbol(data, 16, 8, 3)
% 
%
% Input:
%   data              is the raw time series. 
%   N                 is the length of sliding window (use the length of the raw time series
%                     instead if you don't want to have sliding windows)
%   n                 is the number of symbols in the low dimensional approximation of the sub sequence.
%   alphabet_size     is the number of discrete symbols. 2 <= alphabet_size <= 20, although 
%                     alphabet_size = 2 is a special "useless" case.
%
% Output:
%   symbolic_data:    matrix of symbolic data (no-repetition).  If consecutive subsequences
%                     have the same string, then only the first occurrence is recorded, with
%                     a pointer to its location stored in "pointers"
%   pointers:         location of the first occurrences of the strings
%
% The variable "win_size" is assigned to N/n, this is the number of data points on the raw 
% time series that will be mapped to a single symbol, and can be imagined as the 
% "compression rate".
%
% The symbolic data is returned in "symbolic_data", with pointers to the subsequences  
%
%
% 
%
% Copyright (c) 2003, Eamonn Keogh, Jessica Lin, Stefano Lonardi, Pranav Patel, Li Wei.  All rights reserved.


function [symbolic_data, pointers] =  timeseries2symbol(data, N, n, cutlines, normalize, NR_opt)

    if nargin < 4
        disp('usage: timeseries2symbol(data, window_len, num_segment, cutlines, [normalize_data_option], [numerosity_reduction_option]');
        return;
    end

    norm_thresh = 0.001;

    if nargin < 6
        NR_opt = 2;
    end
    if nargin < 5
        normalize = 1;
    end

    win_size = floor(N/n);                      % win_size is the number of data points on the raw time series that will be mapped to a single symbol
    data_len = max(size(data));

    pointers      = [];                         % Initialize pointers
    sliding_win_num = data_len - (N -1);
    symbolic_data = zeros(data_len - (N -1),n);	% Initialize symbolic_data with a void string, it will be removed later
    all_string    = zeros(data_len-N+1,n);	%#ok<*NASGU>

    % Scan accross the time series, extract subsequences, and convert them to strings
    j = 0;
    for i = 1 : data_len - (N -1)

        % Remove the current subsection
        sub_section = data(i:i + N -1); 

        % Z-normalize it
        if normalize == 1
            sub_section_mean = mean(sub_section);
            sub_section_std = std(sub_section);
            if sub_section_std>norm_thresh
                sub_section = (sub_section - sub_section_mean)/sub_section_std;
            else
                sub_section = sub_section - sub_section_mean;
            end
        end

        PAA = tsPAA(sub_section,n);

        current_string = map_to_string(PAA,cutlines); % Convert the PAA to a string

        % No numerosity reduction: record everything
        if NR_opt == 1
            j = j+1;
            symbolic_data(j,:)	= current_string;     %#ok<*AGROW> % ... add it to the set...
            pointers(j)	= i;                                   % ... and add a new pointer

        % With numerosity reduction: record a string only if it differs from its leftmost neighbor
        elseif NR_opt == 2        
            if ~all(current_string == symbolic_data(end,:))             % If the string differs from its leftmost neighbor...
                j = j+1;
                symbolic_data(j,:)	= current_string;     % ... add it to the set...
                pointers(j)        = i;                      % ... and add a new pointer
            end

        % Advanced numerosity reduction: record a string only if its mindist to the last recorded
        % string > 0
        elseif NR_opt == 3
            % always record the first string
            if i == 1
                j = j+1;
                symbolic_data(j,:)	=  current_string;     % ... add it to the set...
                pointers(j)         = i;                      % ... and add a new pointer

            % subsequent strings
            else            
                % we only need to check if two sliding windows have different strings (if they are
                % the same then their mindist is 0)
                if ~all(current_string == symbolic_data(end,:))            % If the string differs from its leftmost neighbor...                    
                    % Here we're doing a simplified version of mindist. Since we're only interested
                    % in knowing if the distance of two strings is 0, we can do so without any extra
                    % computation. Since only adjacent symbols have distance 0, all we have to
                    % do is check if any two symbols are non-adjacent
                    if any(abs(symbolic_data(end,:) - current_string) > 1)
                        j = j+1;
                        symbolic_data(j,:)	=  current_string;     % ... add it to the set...
                        pointers(j)         = i;                      % ... and add a new pointer
                    end
                end
            end

        else


            % we only need to check if two sliding windows have different strings (if they are
            % the same then their mindist is 0)
            if ~all(current_string == symbolic_data(end,:))            % If the string differs from its leftmost neighbor...        
                if any(abs(symbolic_data(end,:) - current_string) > 1)
                    if ~all(sign(diff(current_string)) >= 0) && ~all(sign(diff(current_string)) <= 0)
                        j = j+1;
                        symbolic_data(j,:)	= current_string;     % ... add it to the set...
                        pointers(j)         = i;                      % ... and add a new pointer
                    end
                end
            end
        end

    end

    % Delete the first element, it was just used to initialize the data structure
    symbolic_data(j+1:end,:) = [];    
end

%%

function string = map_to_string(PAA,cutlines)

    string = zeros(1,length(PAA));

    cut_points = [-Inf cutlines];

    for i = 1 : length(PAA)    
        string(i) = sum( (cut_points <= PAA(i)), 2 );         % order is now: a = 1, b = 2, c = 3..
    end

end
