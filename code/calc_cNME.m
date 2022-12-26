input_data = importdata(filename);
var_names = strip(input_data.textdata,'"');
input_data = input_data.data;
out_mat = Get_cDMI(input_data,var_names,var_names,nz,p,k);
out = array2table(out_mat,'VariableNames',var_names,'RowNames',var_names);
writetable(out,[output,'.csv'],'WriteRowNames',true);