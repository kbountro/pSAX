function y = d2b(x,c)

% Convert a decimanl number into a binary array
% 
% Similar to dec2bin but yields a numerical array instead of a string and is found to
% be rather faster

c = floor(log2(max(x)))+1; % Number of divisions necessary ( rounding up the log2(x) )
y = zeros(length(x),c); % Initialize output array
for i = 1:c
    r = floor(x / 2);
    y(:,c+1-i) = x - 2*r;
    x = r;
end
