function [data,pure_data] = generate_time_series(n,d,pdf,type,sigma,normalize)

    if nargin < 6
        normalize = 1;
    end

    if strcmp(type,'walk')
        data = random_walk(n,d,pdf,sigma);
        pure_data = [];   % placeholder
        
    elseif strcmp(type,'wave')
        f1 = 25; f2 = 66; f3 = 100;
        data = sin(2*pi*f1*(1:n)/8000)+sin(2*pi*f2*(1:n)/8000);%+sin(2*pi*f3*(1:n)/1000);
        pure_data = data;
        if any(pdf==["normal","lognormal"])
            data = data + random(pdf,0,sigma,1,n);
        elseif any(pdf==["poisson","exp"])
            data = data + ( 2*(rand(1,n)>0.5)-1 ).*random(pdf,1,1,n);
        elseif pdf=="laplace"
            data = data + laprnd(1,n,0,sigma);
        end 

    end

    % Z-normalization
    if normalize == 1
        data = ((data' - mean(data'))./std(data'))';
    end
    
end


%------------------------------------------------------------------------------------------
% Make random walk data
%------------------------------------------------------------------------------------------

function r = random_walk(n,d,pdf,sigma)
    % r = random_walk(n)
    % n: length of random walk time series
    % d: dimensionality of time series data
    % 
    % This is the continuous analog of symmetric random walk, each increment y(s+t)-y(s) is 
    % Gaussian with distribution N(0,t^2) and increments over disjoint intervals are independent. 
    % It is typically simulated as an approximating random walk in discrete time. 
    
    for i=1:d
        
        if any(pdf==["normal","lognormal"])
            sigma = 1;
            step_gen = str2func(['@(mu,sigma,M,N)' 'random(''' pdf ''',mu,sigma,M,N)']);
            steps = ( 2*(rand(1,n-1)>0.5)-1 ).*step_gen(0,sigma,1,n-1);
        elseif any(pdf==["poisson","exp"])
            step_gen = str2func(['@(mu,M,N)' 'random(''' pdf ''',mu,M,N)']);
            steps = ( 2*(rand(1,n-1)>0.5)-1 ).*step_gen(1,1,n-1);
        elseif pdf=="laplace"
            steps = laprnd(1,n-1,0,sigma);
        end
        
        r(i,:) = [0 cumsum(steps)];
        
    end
    
    r = r./sqrt(d);
end