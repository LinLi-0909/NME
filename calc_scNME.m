GEM = importdata(filename);
cellnames = GEM.textdata(2:end,1);
gene_names = GEM.textdata(1,2:end);
gene_names = strip(gene_names,'"');
input_data = GEM.data;
if(species == 'human')
    human_tf = importdata('../Data/human_tf_trrust.txt');
    regulators = intersect(gene_names,human_tf);
else
    mouse_tf = importdata('../Data/mouse_tf_trrust.txt');
    regulators = intersect(gene_names,mouse_tf);
end
out_mat = Get_cDMIpar(input_data,regulators,gene_names,nz,p,k,nworkers);
out = array2table(out_mat,'VariableNames',var_names,'RowNames',var_names);

[nme,bnet,regulon] = Get_entropy(GEM,net,regulators,gene_names,p_cutoff);
tnet = array2table(nme,'VariableNames',regulon,'RowNames',cellnames);
writetable(tnet,[output,'_NFmatrix.csv'],'WriteRowNames',true)
nbnet = array2table(bnet,'VariableNames',gene_names,'RowNames',regulon);
writetable(nbnet,[output,'_binarynet.csv'],'WriteRowNames',true)

