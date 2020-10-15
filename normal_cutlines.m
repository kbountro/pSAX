function cutlines = normal_cutlines(alphabet_size)
% Compute equiprobable intervals under the standard Gaussian curve
%
% Author: Konstantinos Bountrogiannis
% Contact: kbountrogiannis@gmail.com
% Date: July 2020

    x = 1/alphabet_size:1/alphabet_size:1-1/alphabet_size;
    cutlines = norminv(x,0,1);
    
end













