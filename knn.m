function [dis] = knn(x,p)
x = x(:);
[dis,~] = knnsearch(x,x,'K', p+1,'Distance','euclidean');
dis = dis(:,2:p+1);
end

