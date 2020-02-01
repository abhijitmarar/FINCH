clc;
clear all;
close all;

wave = 0.0005; %wavelength
f_o = 3.000;   %focal length of objective
D_SLM = 7.68;   %diameter of SLM
f_d = 750;     %focal length of difractive spherical lens (fd = 2*zh is optimal)
d1 = 450;      %distance between objective and SLM
z_h = 600;     %distance between SLM and CCD aperture
NA_in = (7.68/6); 
delta_min = zeros(1,11);
i = 1;


for z_s = 2.995:.0001:3.005
    if (abs(z_s) == 3.000)
        %Simplification f_o = z_s
        z_r = abs(z_h - f_d);
        trans_mag = z_h/f_o;
        delta_min(i) = z_r/(trans_mag*D_SLM);
        
    else
        f_e = (z_s*f_o)/(f_o-z_s);
        f_1 = (f_d*(f_e+d1))/(f_d-(f_e+d1));
        z_r = abs(((f_1+z_h)*(f_e+d1+z_h))/(f_1-f_e-d1));
        trans_mag = (z_h*f_e)/(z_s*(f_e+d1));
        delta_min(i) = z_r/(trans_mag*D_SLM);
       
    end
    i = i+1;
end

z_s = 2.995:.0001:3.005;
plot(z_s,delta_min);