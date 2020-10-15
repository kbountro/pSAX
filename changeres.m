function data_out = changeres(data,res)
% Change input's length (stretch or shrink)
%
% Input:    data     : input data
%           res      : output's length
% Output:   data_out : re-sized data

N = length(data);


    temp = zeros(res, N);
    for j = 1 : res
        temp(j, :) = data;
    end
    expanded_data = reshape(temp, 1, N*res);
    data_out = mean(reshape(expanded_data, N, res));

end

