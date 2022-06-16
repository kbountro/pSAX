% Descritization of the PAA sequence with custom cutlines
%
% Author: Konstantinos Bountrogiannis
% Contact: kbountrogiannis@gmail.com
% Date: June 2022

function string = map_to_string(PAA,cutlines)

    string = zeros(1,length(PAA));

    cut_points = [-Inf cutlines];

    for i = 1 : length(PAA)    
        string(i) = sum( (cut_points <= PAA(i)), 2 );         % order is now: a = 1, b = 2, c = 3..
    end

end