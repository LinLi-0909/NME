%%
% case1_low
clear
N = 2000;
p=15;
beta_xz = 0.5;beta_zy = 0.5;
r_xz = 0.01;r_zy = 0.01;
e1 =  normrnd(0,1, N,1); 
e2 = normrnd(0,1, N,1);
x =  normrnd(0,1, N,1);
z = beta_xz.*abs(x) -0.5+ r_xz.*e1;
y = beta_zy.*abs(2*z) + r_zy.*e2;
case1(:,1) = x;
case1(:,2) = y;
case1(:,3) = z;
out = zeros(3,3);
out(1,2)= cDMI(x,y, z, p, 3);
out(2,1)= cDMI(y,x, z, p, 3);
out(3,2)= cDMI(z,y, x, p, 3);
out(2,3)= cDMI(y,z, x, p, 3);
out(1,3)= cDMI(x,z, y, p, 3);
out(3,1)= cDMI(z,x, y, p, 3);
csvwrite('case1_low.csv',case1)
%%
% case2_low
clear
N = 2000;
beta_xz = 0.5;beta_zy = 0.5; 
r_xyz =0.01;
p=15;
e = normrnd(0,1, N,1);  
x = normrnd(0,1, N,1);
y = normrnd(0,1, N,1);
z = beta_xz.*abs(x).^(1/2) +beta_zy*abs(y).^(1/2) + r_xyz.*e;
case2(:,1) = x;
case2(:,2) = y;
case2(:,3) = z;
out = zeros(3,3);
out(1,2)= cDMI(x,y, z, p, 3);
out(2,1)= cDMI(y,x, z, p, 3);
out(3,2)= cDMI(z,y, x, p, 3);
out(2,3)= cDMI(y,z, x, p, 3);
out(1,3)= cDMI(x,z, y, p, 3);
out(3,1)= cDMI(z,x, y, p, 3);
csvwrite('case2_low.csv',case2)

%%
% case3_low
clear
N = 2000;
beta_xz = 0.5;beta_yz = 0.5; 
r_xz =0.01; r_yz=0.01; p=15;
e1 = normrnd(0,1, N,1);  
e2 = normrnd(0,1, N,1);
z =  normrnd(0,1, N,1);
x = beta_xz.*abs(z+0.5) + r_xz.*e1;
y = beta_yz*abs(z -0.5)+r_yz.*e2;
case3(:,1) = x;
case3(:,2) = y;
case3(:,3) = z;
out = zeros(3,3);
out(1,2)= cDMI(x,y, z, p, 3);
out(2,1)= cDMI(y,x, z, p, 3);
out(3,2)= cDMI(z,y, x, p, 3);
out(2,3)= cDMI(y,z, x, p, 3);
out(1,3)= cDMI(x,z, y, p, 3);
out(3,1)= cDMI(z,x, y, p, 3);
csvwrite('case3_low.csv',case3)

%%
% case4_low
clear
N = 2000;
beta_xz = 0.5;beta_yz = 0.5; beta_xy = 0.5;
r_xz =0.01; r_yz=0.01; p=15;
e1 = normrnd(0,1, N,1);  
e2 = normrnd(0,1, N,1);
z =  normrnd(0,1, N,1);
x = beta_xz.*abs(z)-0.5+r_xz.*e1;
y = beta_yz*(z.^2)+r_yz.*e2 + beta_xy*(x-1).^2;
case4(:,1) = x;
case4(:,2) = y;
case4(:,3) = z;
out = zeros(3,3);
out(1,2)= cDMI(x,y, z, p, 3);
out(2,1)= cDMI(y,x, z, p, 3);
out(3,2)= cDMI(z,y, x, p, 3);
out(2,3)= cDMI(y,z, x, p, 3);
out(1,3)= cDMI(x,z, y, p, 3);
out(3,1)= cDMI(z,x, y, p, 3);
csvwrite('case4_low.csv',case4)

%%
% case5_low
clear
N = 2000;
p=15;
beta_xz = 0.5;beta_zy = 0.5;beta_yx = 0.5;r_xyz = 0.01;
x = zeros(2000,1);
y = zeros(2000,1);
z = zeros(2000,1);
e1 = normrnd(0,1,N,1);
e2 = normrnd(0,1,N,1);
e3 = normrnd(0,1,N,1);
x(1,1) = normrnd(0,1,1,1);
y(1,1) = normrnd(0,1,1,1);
z(1,1) = normrnd(0,1,1,1);
for i = 2:2000
z(i,1) = beta_xz.*abs(2*x(i-1,1))-0.5+ r_xyz.*e2(i,1);
y(i,1) = beta_xz.*abs(2*z(i-1,1))+ r_xyz.*e3(i,1);
x(i,1) = beta_xz.*abs(2*y(i-1,1)-1)+ r_xyz.*e1(i,1);
end  
case5(:,1) = x;
case5(:,2) = y;
case5(:,3) = z;
out = zeros(3,3);
out(1,2)= cDMI(x,y, z, p, 3);
out(2,1)= cDMI(y,x, z, p, 3);
out(3,2)= cDMI(z,y, x, p, 3);
out(2,3)= cDMI(y,z, x, p, 3);
out(1,3)= cDMI(x,z, y, p, 3);
out(3,1)= cDMI(z,x, y, p, 3);
csvwrite('case5_low.csv',case5)
%%
% case1_median
clear
N = 2000;
p=15;
beta_xz = 0.5;beta_zy = 0.5;
r_xz = 0.05;r_zy = 0.05;
e1 =  normrnd(0,1, N,1); 
e2 = normrnd(0,1, N,1);
x =  normrnd(0,1, N,1);
z = beta_xz.*abs(x) -0.5+ r_xz.*e1;
y = beta_zy.*abs(2*z) + r_zy.*e2;
case1(:,1) = x;
case1(:,2) = y;
case1(:,3) = z;
out = zeros(3,3);
out(1,2)= cDMI(x,y, z, p, 3);
out(2,1)= cDMI(y,x, z, p, 3);
out(3,2)= cDMI(z,y, x, p, 3);
out(2,3)= cDMI(y,z, x, p, 3);
out(1,3)= cDMI(x,z, y, p, 3);
out(3,1)= cDMI(z,x, y, p, 3);
csvwrite('case1_median.csv',case1)
%%
% case2_median
clear
N = 2000;
beta_xz = 0.5;beta_zy = 0.5; 
r_xyz =0.05;
p=15;
e = normrnd(0,1, N,1);  
x = normrnd(0,1, N,1);
y = normrnd(0,1, N,1);
z = beta_xz.*abs(x).^(1/2) +beta_zy*abs(y).^(1/2) + r_xyz.*e;
case2(:,1) = x;
case2(:,2) = y;
case2(:,3) = z;
out = zeros(3,3);
out(1,2)= cDMI(x,y, z, p, 3);
out(2,1)= cDMI(y,x, z, p, 3);
out(3,2)= cDMI(z,y, x, p, 3);
out(2,3)= cDMI(y,z, x, p, 3);
out(1,3)= cDMI(x,z, y, p, 3);
out(3,1)= cDMI(z,x, y, p, 3);
csvwrite('case2_median.csv',case2)

%%
% case3_median
clear
N = 2000;
beta_xz = 0.5;beta_yz = 0.5; 
r_xz =0.05; r_yz=0.05; p=15;
e1 = normrnd(0,1, N,1);  
e2 = normrnd(0,1, N,1);
z =  normrnd(0,1, N,1);
x = beta_xz.*abs(z+0.5) + r_xz.*e1;
y = beta_yz*abs(z -0.5)+r_yz.*e2;
case3(:,1) = x;
case3(:,2) = y;
case3(:,3) = z;
out = zeros(3,3);
out(1,2)= cDMI(x,y, z, p, 3);
out(2,1)= cDMI(y,x, z, p, 3);
out(3,2)= cDMI(z,y, x, p, 3);
out(2,3)= cDMI(y,z, x, p, 3);
out(1,3)= cDMI(x,z, y, p, 3);
out(3,1)= cDMI(z,x, y, p, 3);
csvwrite('case3_median.csv',case3)

%%
% case4_median
clear
N = 2000;
beta_xz = 0.5;beta_yz = 0.5; beta_xy = 0.5;
r_xz =0.05; r_yz=0.05; p=15;
e1 = normrnd(0,1, N,1);  
e2 = normrnd(0,1, N,1);
z =  normrnd(0,1, N,1);
x = beta_xz.*abs(z)-0.5+r_xz.*e1;
y = beta_yz*(z.^2)+r_yz.*e2 + beta_xy*(x-1).^2;
case4(:,1) = x;
case4(:,2) = y;
case4(:,3) = z;
out = zeros(3,3);
out(1,2)= cDMI(x,y, z, p, 3);
out(2,1)= cDMI(y,x, z, p, 3);
out(3,2)= cDMI(z,y, x, p, 3);
out(2,3)= cDMI(y,z, x, p, 3);
out(1,3)= cDMI(x,z, y, p, 3);
out(3,1)= cDMI(z,x, y, p, 3);
csvwrite('case4_median.csv',case4)

%%
% case5_median
clear
N = 2000;
p=15;
beta_xz = 0.5;beta_zy = 0.5;beta_yx = 0.5;r_xyz = 0.05;
x = zeros(2000,1);
y = zeros(2000,1);
z = zeros(2000,1);
e1 = normrnd(0,1,N,1);
e2 = normrnd(0,1,N,1);
e3 = normrnd(0,1,N,1);
x(1,1) = normrnd(0,1,1,1);
y(1,1) = normrnd(0,1,1,1);
z(1,1) = normrnd(0,1,1,1);
for i = 2:2000
z(i,1) = beta_xz.*abs(2*x(i-1,1))-0.5+ r_xyz.*e2(i,1);
y(i,1) = beta_zy.*abs(2*z(i-1,1))+ r_xyz.*e3(i,1);
x(i,1) = beta_yx.*abs(2*y(i-1,1)-1)+ r_xyz.*e1(i,1);
end  
case5(:,1) = x;
case5(:,2) = y;
case5(:,3) = z;
out = zeros(3,3);
out(1,2)= cDMI(x,y, z, p, 3);
out(2,1)= cDMI(y,x, z, p, 3);
out(3,2)= cDMI(z,y, x, p, 3);
out(2,3)= cDMI(y,z, x, p, 3);
out(1,3)= cDMI(x,z, y, p, 3);
out(3,1)= cDMI(z,x, y, p, 3);
csvwrite('case5_median.csv',case5)
%%
% case1_high
clear
N = 2000;
p=15;
beta_xz = 0.5;beta_zy = 0.5;
r_xz = 0.1;r_zy = 0.1;
e1 =  normrnd(0,1, N,1); 
e2 = normrnd(0,1, N,1);
x =  normrnd(0,1, N,1);
z = beta_xz.*abs(x) -0.5+ r_xz.*e1;
y = beta_zy.*abs(2*z) + r_zy.*e2;
case1(:,1) = x;
case1(:,2) = y;
case1(:,3) = z;
out = zeros(3,3);
out(1,2)= cDMI(x,y, z, p, 3);
out(2,1)= cDMI(y,x, z, p, 3);
out(3,2)= cDMI(z,y, x, p, 3);
out(2,3)= cDMI(y,z, x, p, 3);
out(1,3)= cDMI(x,z, y, p, 3);
out(3,1)= cDMI(z,x, y, p, 3);
csvwrite('case1_high.csv',case1)
%%
% case2_high
clear
N = 2000;
beta_xz = 0.5;beta_zy = 0.5; 
r_xyz =0.1;
p=15;
e = normrnd(0,1, N,1);  
x = normrnd(0,1, N,1);
y = normrnd(0,1, N,1);
z = beta_xz.*abs(x).^(1/2) +beta_zy*abs(y).^(1/2) + r_xyz.*e;
case2(:,1) = x;
case2(:,2) = y;
case2(:,3) = z;
out = zeros(3,3);
out(1,2)= cDMI(x,y, z, p, 3);
out(2,1)= cDMI(y,x, z, p, 3);
out(3,2)= cDMI(z,y, x, p, 3);
out(2,3)= cDMI(y,z, x, p, 3);
out(1,3)= cDMI(x,z, y, p, 3);
out(3,1)= cDMI(z,x, y, p, 3);
csvwrite('case2_high.csv',case2)

%%
% case3_high
clear
N = 2000;
beta_xz = 0.5;beta_yz = 0.5; 
r_xz =0.01; r_yz=0.01; p=15;
e1 = normrnd(0,1, N,1);  
e2 = normrnd(0,1, N,1);
z =  normrnd(0,1, N,1);
x = beta_xz.*abs(z+0.5) + r_xz.*e1;
y = beta_yz*abs(z -0.5)+r_yz.*e2;
case3(:,1) = x;
case3(:,2) = y;
case3(:,3) = z;
out = zeros(3,3);
out(1,2)= cDMI(x,y, z, p, 3);
out(2,1)= cDMI(y,x, z, p, 3);
out(3,2)= cDMI(z,y, x, p, 3);
out(2,3)= cDMI(y,z, x, p, 3);
out(1,3)= cDMI(x,z, y, p, 3);
out(3,1)= cDMI(z,x, y, p, 3);
csvwrite('case3_high.csv',case3)

%%
% case4_high
clear
N = 2000;
beta_xz = 0.5;beta_yz = 0.5; beta_xy = 0.5;
r_xz =0.1; r_yz=0.1; p=15;
e1 = normrnd(0,1, N,1);  
e2 = normrnd(0,1, N,1);
z =  normrnd(0,1, N,1);
x = beta_xz.*abs(z)-0.5+r_xz.*e1;
y = beta_yz*(z.^2)+r_yz.*e2 + beta_xy*(x-1).^2;
case4(:,1) = x;
case4(:,2) = y;
case4(:,3) = z;
out = zeros(3,3);
out(1,2)= cDMI(x,y, z, p, 3);
out(2,1)= cDMI(y,x, z, p, 3);
out(3,2)= cDMI(z,y, x, p, 3);
out(2,3)= cDMI(y,z, x, p, 3);
out(1,3)= cDMI(x,z, y, p, 3);
out(3,1)= cDMI(z,x, y, p, 3);
csvwrite('case4_high.csv',case4)

%%
% case5_high
clear
N = 2000;
p=15;
beta_xz = 0.5;beta_zy = 0.5;beta_yx = 0.5;r_xyz = 0.1;
x = zeros(2000,1);
y = zeros(2000,1);
z = zeros(2000,1);
e1 = normrnd(0,1,N,1);
e2 = normrnd(0,1,N,1);
e3 = normrnd(0,1,N,1);
x(1,1) = normrnd(0,1,1,1);
y(1,1) = normrnd(0,1,1,1);
z(1,1) = normrnd(0,1,1,1);
for i = 2:2000
z(i,1) = beta_xz.*abs(2*x(i-1,1))-0.5+ r_xyz.*e2(i,1);
y(i,1) = beta_xz.*abs(2*z(i-1,1))+ r_xyz.*e3(i,1);
x(i,1) = beta_xz.*abs(2*y(i-1,1)-1) + r_xyz.*e1(i,1);
end  
case5(:,1) = x;
case5(:,2) = y;
case5(:,3) = z;
out = zeros(3,3);
out(1,2)= cDMI(x,y, z, p, 3);
out(2,1)= cDMI(y,x, z, p, 3);
out(3,2)= cDMI(z,y, x, p, 3);
out(2,3)= cDMI(y,z, x, p, 3);
out(1,3)= cDMI(x,z, y, p, 3);
out(3,1)= cDMI(z,x, y, p, 3);
csvwrite('case5_high.csv',case5)