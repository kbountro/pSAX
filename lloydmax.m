function [C,b] = lloydmax(f,x,ncodewords,data,init,lock_bounds)
% Quantization using Lloyd-Max quantizer.
% Inputs:   f(x): density function
%           init: Initialization. Either initial bounds or initial
%           codewords. Can be empty.
%           lock_bounds: Locks initial bounds and optimizes only codewords.
%           Effectively, it computes the centroids of the given bounds.
%           
% Output:   C: codewords
%           b: boundaries that define the partition
%
% Author: Konstantinos Bountrogiannis
% Contact: kbountrogiannis@gmail.com
% Date: July 2020

    if nargin<6 || isempty(lock_bounds)
        lock_bounds = 0;
    end
    
    %% Initialization
    
    d = min(size(x));
    
    % Fix dimensions of 'f' and 'x'
    if d>1
        f = f(:)';
        for k = 1:d
            % first clear x to reset its dimensions
            if k==1, x_temp = x; clear x; end
            x(k,:) = x_temp(:,k);
        end
        clear x_temp;
    end
    
    % Determine 'initial' input
    if ~isempty(init)
        if length(init) == ncodewords
            init_bounds = 0;
            init_codebook = 1;
            C = init;
        elseif length(init) == ncodewords-1
            init_bounds = 1;
            init_codebook = 0;
            b = init;
        elseif length(init)<ncodewords % Check if initial codebook does not match the desired size
            % and fill-in if needed
            [~,idx] = min( abs(init-mean(init)) );
            v = init(idx);
            space = ncodewords-length(init);
            C = [init(1:idx) ones(1,space)*v init(idx+1:end)];
            init_bounds = 0;
            init_codebook = 1;
        end
    else
        init_bounds = 0;
        init_codebook = 0;
    end
    
    % Initialize C and b
    if ~init_codebook && ~isempty(data)
        % k-means++ initialization.
        [~,C] = kmeanspp(data,ncodewords);
    end
    prev_C = zeros(1,ncodewords);
    prev_b = zeros(1,ncodewords-1);
    
    %% Main algorithm
    
    stop = 0;
    iter = 0;
    while ~stop && iter<100
        if (iter>1 || ~init_bounds) && ~lock_bounds
            b = (C(1:1:end-1)+C(2:1:end))/2;
        end
        for i = 1:ncodewords
            switch i
                case 1
                    b1 = 1;
                    [~,b2] = min(abs(b(i)-x));
                    if b2 == b1, b2 = b1+1; end % Bounds must not coincide
                case ncodewords
                    [~,b1] = min(abs(b(i-1)-x));
                    b2 = length(x);
                    if b2 == b1, b1 = b2-1; end % Bounds must not coincide
                otherwise
                    [~,b1] = min(abs(b(i-1)-x));
                    [~,b2] = min(abs(b(i)-x));
                    if b2 == b1 % Bounds must not coincide
                        if b1>1
                            b1 = b2-1;
                        elseif b2<length(x)
                            b2 = b1+1;
                        end
                    end
            end
    
            I1 = trapz(x(b1:b2),x(b1:b2).*f(b1:b2),2);
            I2 = trapz(x(b1:b2),f(b1:b2),2);
            C(i) = I1/I2;
            if ( isnan(C(i)) || isinf(C(i)) )
                f(b1:b2) = f(b1:b2)+realmin; f = f/trapz(x,f,2);
                I1 = trapz(x(b1:b2),x(b1:b2).*f(b1:b2),2);
                I2 = trapz(x(b1:b2),f(b1:b2),2);
                C(i) = I1/I2;
                if isnan(C(i)) || isinf(C(i))
                    error('Error: lloydmax returned NaN or Inf codewords');
                end
            end
        end

        stop = all(prev_C==C) && all(prev_b==b);
        
        prev_b = b;
        prev_C = C;
        
        iter = iter+1;
    end
    
end

