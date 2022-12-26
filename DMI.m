function [c] = DMI(x, y, p, k)
%DMI: compute the direct entropy causality from x to y.
% x: N samples
% y: N samples
% p: the number of neighbors of x 
% k: the k-th nearest number to use in calculating entropy, at least 2.
%Use the kNN method to estimate the MI(YNN, Y(xNN)), where the dimension of YNN is p.
%kNN formular MI(YNN, Y(xNN)) = psi(k) - <psi(nYNN+1) + psi(nY(xNN)+1)> + psi(N+1).
N= size(x,1); 
% define the neighbors of x 
% xi should be excluded from the neighbors of xi
dx = pdist2(x, x,'euclidean');
dy = pdist2(y,y,'euclidean');
dx(logical(eye(size(dx))))=exp(30);
dy(logical(eye(size(dy))))=exp(30);
[~,idx]=sort(dx,2);
[~,idy]=sort(dy,2);
XNNid= idx(:,1:p);
YNNid= idy(:,1:p);
YxNN = y(XNNid);  % N*p 
YNN = y(YNNid);  % N*p
%Calculate mutual information of prediction
% kNN in (X, YpNN) space
 nYNN = zeros(N, 1); nYxNN = zeros(N, 1);
 [~, D] = knnsearch([YNN, YxNN], [YNN,YxNN], 'K', k+1, 'Distance', 'chebychev'); %kth NN
 halfepsilon = D(:,k+1); %the distance of kNN = epsilon/2
 nYNN = sum((pdist2(YNN, YNN,'chebychev')<halfepsilon),2)-1; %YNN exclude the point xi
 nYxNN= sum((pdist2(YxNN, YxNN, 'chebychev')<halfepsilon),2)-1;
 idN = (nYNN~=-1)&(nYxNN~=-1);
 c = abs(psi(k) - mean(psi(nYNN(idN)+1)) - mean(psi(nYxNN(idN)+1)) + mean(psi(N+1)));

   [~, Da] = knnsearch([YNN, YNN], [YNN,YNN], 'K', k+1, 'Distance', 'chebychev'); %kth NN
 halfepsilon_a = Da(:,k+1); %the distance of kNN = epsilon/2
  nYNN_a = sum((pdist2(YNN, YNN, 'chebychev')<halfepsilon_a),2)-1; %YNN exclude the point xi
  a = psi(k) - 2*mean(psi(nYNN_a+1)) + mean(psi(N+1));
  e=0.00001;
 out = c/a; 
end
