clear all;
clc;

z_s = 2.995;
f_o = 3.000;   %focal length of objective
%d1 = 0;      %distance between objective and SLM
z_h = 300;     %distance between SLM and CCD aperture
mag = zeros(1,31);
counter = 1;
for d1 = 0:10:300
    f_e = (z_s*f_o)/(f_o-z_s);
    mag(counter) = (z_h*f_e)/(z_s*(f_e+d1));
    counter = counter+1;
end
%z_s = 2.995:0.0001:3.005;
d1 = 0:10:300;
plot(d1,mag)