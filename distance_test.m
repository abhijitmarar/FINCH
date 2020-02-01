clc;
clear all;
d2 = 4e+5;                  %SLM to CCD
z = 10;                     %defocus
d1 = 1e+5:1e+3:4e+5;        %obj to SLM
%a = 8e+5;                   %f_doe
f = 9e+3;                   %focal length
check = 2.995:2.0e-6:3.005; 
y = zeros(1,301);
for a=2.995:0.002:3.005
    num = -400*((600*a)-1801)*((600*a-1799));
    den = 3;
    z_r = num./den;
end