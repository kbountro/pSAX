function [C,b] = lloydmax(f,x,ncodewords,data,init,lock_bounds)
% Quantization using Lloyd-Max quantizer.
% Inputs:   f: pdf
%           x: data
%           init: Initialization. Either initial bounds or initial codewords. Can be empty.
%           lock_bounds: (0 or 1). Locks initial bounds and optimizes only codewords.
%                        Efffectively, it computes the centroids of the given bounds. 
%                        Normally, this is set to 0.
%           
% Output:   C: Codeword values
%           b: boundaries of quantization intervals
%
% Author: Konstantinos Bountrogiannis
% Contact: kbountrogiannis@gmail.com
% Date: June 2022

    if nargin<6
        lock_bounds = 0;
    end
    
    %% Initialization
    tol = 1e-7;
    tol2 = eps*(max(x)-min(x));
    
    % Determine 'initial' input
    C = [];
    if nargin>4
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
        C = sort(C);
    end
    prev_C = zeros(1,ncodewords);
    prev_b = zeros(1,ncodewords-1);
    
    %% Main algorithm
    
    stop = 0;
    iter = 0;
    if isempty(C), C = mod(x(1)+(x(end)-x(1))*rand(1,2),x(end)); end
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
                    if b2 <= b1 % Bounds must not be in reverse order
                        if b1>1
                            b1 = b2-1;
                        elseif b2<length(x)
                            b2 = b1+1;
                        end
                    end
            end
    
            b2 = min(length(x),b2); % this activates only when length(x)==1
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

%         stop = sum(all(prev_C==C)) && sum(all(prev_b==b));
        move = mean(abs(prev_C-C));
        if move > tol2
            rel_dist = move/mean(abs(C));
        else
            rel_dist = move;
        end
        stop = rel_dist<tol && rel_dist<tol2;
        
        prev_b = b;
        prev_C = C;
        
        iter = iter+1;
    end
    
end

