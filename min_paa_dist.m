% THIS IS A MINOR MODIFICATION OF THE min_dist.m FILE, WHICH IS
% (C) Eamonn Keogh, Jessica Lin, Stefano Lonardi, Pranav Patel, Li Wei.
% THE PRESENT FILE IS WRITTEN BY Konstantinos Bountrogiannis.

function dist = min_paa_dist(str1_PAA, str2_SAX, alphabet_size, compression_ratio, cutlines)
% This function computes the minimum (lower-bounding) distance between a
% symbolic and a PAA sequence.
% The vectors should have equal length.
%   
%   Input:
%       str1_PAA: PAA sequence
%       str2_SAX: symbolic sequence
%       alphabet_size: alphabet size used to construct the symbolic sequence
%       compression_ratio: original_data_len / symbolic_len
%   Output:
%       dist: lower-bounding distance
%
%   usage: dist = min_dist(str1_PAA, str2_SAX, alphabet_size, compression_ratio, cutlines)
%

    if (length(str1_PAA) ~= length(str2_SAX))
        disp('error: the strings must have equal length!');
        return;
    end
    
%     if (any(str1_PAA > alphabet_size) || any(str2_SAX > alphabet_size))
%         disp('error: some symbol(s) in the string(s) exceed(s) the alphabet size!');
%         return;
%     end
    
    dist = 0; %#ok<*NASGU>

    str_len = length(str1_PAA);
    
    for i = 1 : str_len
        if str2_SAX(i) > 1
            if cutlines(str2_SAX(i)-1) > str1_PAA(i)
                dist = dist + (cutlines(str2_SAX(i)-1) - str1_PAA(i))^2;
                continue;
            end
        end
        if str2_SAX(i) < alphabet_size
            if cutlines(str2_SAX(i)) < str1_PAA(i)
                dist = dist + (cutlines(str2_SAX(i)) - str1_PAA(i))^2;
            end
        end
    end
    
    dist = sqrt(compression_ratio * dist);
