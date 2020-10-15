% k-means++
% Taken from the k-means file of Laurent S.:
% (https://www.mathworks.com/matlabcentral/fileexchange/28804-k-means)

function [L,C] = kmeanspp(X,k)

L = [];
L1 = 0;
% k = min(length(X),k);

iter = 0;
while length(unique(L)) ~= k && iter<10
    
    % The k-means++ initialization.
    C = X(:,1+round(rand*(size(X,2)-1)));
    L = ones(1,size(X,2));
    for i = 2:k
        D = X-C(:,L);
        D = cumsum(sqrt(dot(D,D,1)));
        if D(end) == 0, C(:,i:k) = X(:,ones(1,k-i+1)); return; end
        C(:,i) = X(:,find(rand < D/D(end),1));
        [~,L] = max(2*C'*X-dot(C,C,1).');
    end
%     
%     % The k-means algorithm.
%     while any(L ~= L1)
%         L1 = L;
%         for i = 1:k, l = L==i; C(:,i) = sum(X(:,l),2)/sum(l); end
%         [~,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C,1).'),[],1);
%     end
    iter = iter+1;
end

C = sort(C);

