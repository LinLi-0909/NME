function [out] = Get_cDMI(input_data,regulators,gene_names,nz,p,l)
%Get_cDMI: calculate the weight matrix of GRN graph
%GEM: variables gene expression matrix samples x genes m*n
%regulators: tf gene cell variable included in gene names with nreg *1
%gene_names: all gene names cell variable corresponding with GEM column with ngene*1
% nz: number of z,choose top nz of the third varible
% condition gene contain all genes
% cDMI: compute the direct entropy causality from x to y.

[~,pos] = ismember(regulators,gene_names);
pos(pos==0)=[];

[m,n] = size(input_data);
dc_GEM = mat2cell(input_data,m,ones(1,n));
A_id = cellfun(@(x) knn(x,p),dc_GEM','UniformOutput',false);
% if regulators
B_id = A_id(pos,:);
%define the index of z 
R = corr(input_data);
[~,idr] = sort(sum(abs(R),2),'descend');
cDMI = zeros(length(regulators),n);
for k = 1:nz
  disp(['The conditional gene ',num2str(k)]);
  z = input_data(:,idr(k));
  ZNN = z(A_id{idr(k),1});
  [~, D2] = knnsearch([ZNN,ZNN],[ZNN,ZNN], 'K', l+1, 'Distance', 'chebychev'); %kth NN
  halfepsilon2 = D2(:,l+1); %the distance of kNN = epsilon/2
  nZNN_a = sum((pdist2(ZNN, ZNN, 'chebychev')<halfepsilon2),2)-1; %YNN exclude the point xi
  ZxNN = cellfun(@(x) mapping(x,z),B_id,'UniformOutput',false);
  D1 = cellfun(@(x) knn_estimation(x,ZNN,l),ZxNN,'UniformOutput',false);
  %halfepsilon1 = D1(:,l+1); %the distance of kNN = epsilon/2
  halfepsilon1 = cellfun(@(x) halfepsilon(x,l),D1,'UniformOutput',false);
  nZxNN = cellfun(@(x,y)dist_sum(x,y),ZxNN,halfepsilon1,'UniformOutput',false);
  nZNN= cellfun(@(y)dist_sum(ZNN,y),halfepsilon1,'UniformOutput',false);
  YNN = cellfun(@(x,y) mapping(x,y),A_id,dc_GEM','UniformOutput',false);
  YzNN = cellfun(@(x) mapping(A_id{idr(k),1},x),dc_GEM','UniformOutput',false);
  D3 = cellfun(@(x,y) knn_estimation(x,y,l),YzNN,YNN,'UniformOutput',false);
  halfepsilon3 = cellfun(@(x) halfepsilon(x,l),D3,'UniformOutput',false);
  nYzNN = cellfun(@(x,y)dist_sum(x,y),YzNN,halfepsilon3,'UniformOutput',false);
  nYNN =  cellfun(@(x,y)dist_sum(x,y),YNN,halfepsilon3,'UniformOutput',false);  
 for    j= 1:n 
  disp(['Compute the gene ',num2str(j)]);
  y_exp = input_data(:,j);
  YxNN = cellfun(@(x) mapping(x,y_exp), B_id,'UniformOutput',false);
  
  D5 = cellfun(@(x) knn_estimation(x,YNN{j,:},l),YxNN,'UniformOutput',false);
  halfepsilon5 = cellfun(@(x) halfepsilon(x,l),D5,'UniformOutput',false);
  nYxNN  =  cellfun(@(x,y)dist_sum(x,y),YxNN,halfepsilon5,'UniformOutput',false);  
  nYNN2= cellfun(@(x) dist_sum(YNN{j,:},x),halfepsilon5,'UniformOutput',false);  
   idN = cellfun(@(x1,x2,x3,x4) adjusted_index(x1,x2,x3,x4,nYNN{j,:},nYzNN{j,:},nZNN_a),nZNN,nZxNN,nYNN2,nYxNN,'UniformOutput',false);  
   tmp1 = cellfun(@(x,y) psi_mean(x,y),nZNN,idN);  
   tmp2 = cellfun(@(x,y) psi_mean(x,y),nZxNN,idN);  
   tmp3 = cellfun(@(y) psi_mean(nZNN_a,y),idN);  
   tmp4 = cellfun(@(x) psi_mean(nYNN{j,:},x),idN);  
   tmp5 = cellfun(@(x) psi_mean(nYzNN{j,:},x),idN);  
   tmp6 = cellfun(@(x,y) psi_mean(x,y),nYNN2,idN);  
   tmp7 = cellfun(@(x,y) psi_mean(x,y),nYxNN,idN);  
   tmp8 = cellfun(@(x) psi_sum(x),idN);  
   I1 = abs(psi(l) - tmp1 - tmp2 +tmp8);
   I2 = abs(psi(l) -2*tmp3 + tmp8);
   I3 = abs(psi(l) -  tmp4 -tmp5+ tmp8);
   I5 =  abs(psi(l) -tmp6  -tmp7 + tmp8);
   TMP(:,j) =max(I5-I1./I2.*I3,0);
 end
if all(all(cDMI==0))
    cDMI=TMP;
else
  cDMI=min(TMP,cDMI);
end
 end
out = cDMI;
for i = 1:length(pos)
    out(i,pos(i)) = 0;
end
%out = array2table(out,'VariableNames',gene_names,'RowNames',gene_names(pos));
end
