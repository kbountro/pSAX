function out = mse(x,y)
% Mean Square Error
%
% Author: Konstantinos Bountrogiannis
% Contact: kbountrogiannis@gmail.com
% Date: July 2020
%

% if ~all(size(x)==size(y))
%     disp('Input sizes do not match!'); return;
% end


if length(x)==length(y)
    out = mean( (x-y).^2, 'all' );
else
    % {Suppose that some divides the other}
    out = 0;
    if max(length(x),length(y))==length(x)
        seg_length = length(x)/length(y);
        for i = 1:length(y)
            out = out + sum( (x((i-1)*seg_length+1:i*seg_length)-y(i)).^2 );
        end
        out = out/length(x);
    else
        seg_length = length(y)/length(x);
        for i = 1:length(x)
            out = out + sum( (y((i-1)*seg_length+1:i*seg_length)-x(i)).^2 );
        end
        out = out/length(y);
    end
end

end

