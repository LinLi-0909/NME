function [out] = dist_sum(x,y)
out = sum((pdist2(x,x,'chebychev')<y),2)-1; %YNN exclude the point xi
end