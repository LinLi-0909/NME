
input_data = importdata(filename);


x = input_data(:,1);
y = input_data(:,2);
out_x2y = DMI(x,y,p,k);
out_y2x = DMI(y,x,p,k);
disp(['The normalized causality from V1 to V2 is:',num2str(out_x2y)]);
disp(['The normalized causality from V2 to V1 is:',num2str(out_y2x)]);



