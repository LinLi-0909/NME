
function [out,bnet,regulators] = Get_entropy(GEM,net,regulators,gene_names,p_value)
% GEM, cell*gene
% net, TF*genes m*n
% p_value
 [m,~] = size(net);
 [~,pos] = ismember(regulators,gene_names);
 pos(pos==0)=[];
 for i = 1:m
 cnet = net(i,:);
 cmean = mean(cnet);
 sd = std(cnet);
 sig = norminv( 1-p_value,cmean,sd);
  cnet(cnet<sig)=0;
 cnet(cnet>=sig)=1;
 regulators{i}=[regulators{i},'(',num2str(sum(cnet)),'g)'];
 bnet(i,:)= cnet;
 tmp  = cnet.*GEM;
 tmp_s = sum(tmp,2);
 tmp_s(tmp_s==0) =1;
 tmp1 = tmp./tmp_s;
 tmp1(tmp1==0)=1;
 out(:,i) = -sum(tmp1.*log(tmp1).*GEM(:,pos(i)),2);
 end
end
