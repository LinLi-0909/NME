function [auroc,aupr] = get_auc(weightmat,regulators,gene_names,gold)
% get the auc num of weighmat predicted
% weightmat: the weight of each edges in GRN
% goldfilepath: the true edge of the GRN
m = length(gold.data);
predict = zeros(m,1);
for i = 1:m
    a = gold.textdata{i,1};
    b = gold.textdata{i,2};
    [~,ida] = intersect(regulators,a);
    [~,idb] = intersect(gene_names,b);
    predict(i,:) = weightmat(ida,idb);
end

pos_num = sum(gold.data == 1);
neg_num = sum(gold.data == 0);

[~, index] = sort(predict,'descend');
ground_truth = gold.data(index);
predict = predict(index);

FPR = zeros(m+1,1);
TPR = zeros(m+1,1);
auroc = 0;
FPR(1) = 0; TPR(1) = 0;
j=2;
for i = 1:(m-1)
    if ((i==1 || ground_truth(i)~=ground_truth(i+1)) && predict(i)~=0) 
    TP = sum(ground_truth(1:i)==1);
    FP = sum(ground_truth(1:i)==0);
    FPR(j) = FP/neg_num;
    TPR(j) = TP/pos_num;
    j=j+1;
    end
end
FPR(j) = 1;
TPR(j)=1;
FPR = FPR(1:j);TPR=TPR(1:j);
for i = 2:j
    auroc = auroc + (FPR(i)-FPR(i-1))*(TPR(i)+TPR(i-1))/2;
end


recall = zeros(m+1,1);
precision = zeros(m+1,1);
aupr = 0;
recall(1) = 0; precision(1) = pos_num/m;
tmp = 1;
for i = 1:m
    TP = sum(ground_truth(1:i)==1);
    precision(i+1) = TP/i;
    recall(i+1) = TP/pos_num;
    if precision(i+1) < precision(i)
        aupr = aupr + precision(i)*(recall(i)-recall(tmp))/2;
        tmp = i+1;
    end
    %aupr = aupr + (precision(i)+precision(i-1))*(recall(i-1)-recall(i))/2;
end


end

