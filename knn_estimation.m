function [out] = knn_estimation(x,ZNN,l)
[~,out]=knnsearch([x,ZNN], [x,ZNN], 'K', l+1, 'Distance', 'chebychev'); %kth NN
end