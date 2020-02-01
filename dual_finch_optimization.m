%close all
clear all

% SIDH with 2 spherical waves

f_o = 3;                                % Focal length of objective (mm)    
NA = 1.4;                               % Numerical aperture of objective
R_o = (f_o*NA)*(55/18);                 % Radius of beam at interferometer using 60x oil objective (mm)
wave = 680e-6;                          % Wavelength of light (mm)
delta_c = 16e-3;                        % Pixel size of camera (mm)
f_d1 = 434.7826 ; f_d2 = 588.2353;          % Focal lengths of mirrors
f_tl = 180;                             % Focal length of tube lens
f_1 = 550;                              % Focal length of L1
s_fac = (f_d2-f_d1)/(f_d2+f_d1);        % s-factor
z_s = 2.990:200e-6:3.010;
d1 = 183;                               % Distance between objective and tube lens
d2 = 730;                               % Distance between tube lens and L1
d3 = 550;                               % Distance between L1 and Interferometer
z_h = (2*f_d1*f_d2)/(f_d1+f_d2);        % Distance between interferometer and camera (maximum overlap)
%z_h = 450;
z_r = zeros(1,length(z_s));
trans_mag = zeros(1,length(z_s));
magnified_hologram = zeros(1,length(z_s)); 
f_e = zeros(1,length(z_s));
f_g = zeros(1,length(z_s));
f_h = zeros(1,length(z_s));
f_j = zeros(1,length(z_s));
gamma = zeros(1,length(z_s));
%Condition to sample finest fringe of Hologram
z_h_min = ((4*R_o*delta_c)/wave);  % Min SLM-camera distance (mm)

%Reconstruction Distance
for i = 1:length(z_s)
    if (z_s(i) == f_o)
        z_r(i) = -((z_h-f_d1)*(z_h-f_d2))/(f_d1-f_d2);
    else
        f_e(i) = (f_o*z_s(i))/(f_o-z_s(i));
        f_g(i) = (f_tl*(f_e(i)+d1))/(f_tl-(f_e(i)+d1));
        f_h(i) = (f_1*(f_g(i)+d2))/(f_1-(f_g(i)+d2));
        gamma(i) = f_h(i)+d3;
        zf1 = gamma(i)*z_h - f_d1*(gamma(i)+z_h); zf2 = gamma(i)*z_h - f_d2*(gamma(i)+z_h);
        z_r(i) = -(zf1*zf2)/(gamma(i)^2*(f_d1-f_d2));
    end
end

%Radius of hologram at Interferometer
y_enter = 0;
theta_enter = NA;
%Element B of Jones Matrix
elem_B_term5 = (d1*(1-(d3/f_1))*(1-(d2/f_tl)));
elem_B_term6 = (d3*(1-(d1/f_tl)));
elem_B_term7 = (d2*(1-(d3/f_1)));
elem_B = zeros(1,length(z_s));
y_exit = zeros(1,length(z_s));
for i = 1: length(z_s)
    elem_B_term1 = (1-(d1/f_o))*(1-(d3/f_1))*(1-(d2/f_tl))*(z_s(i)); 
    elem_B_term2 = (d3/f_tl)*(1-(d1/f_o))*(z_s(i));
    elem_B_term3 = (d2/f_o)*(1-(d3/f_1))*(z_s(i));
    elem_B_term4 = (d3/f_o)*(z_s(i));  
    elem_B(i) = elem_B_term1-elem_B_term2-elem_B_term3-elem_B_term4+elem_B_term5+elem_B_term6+elem_B_term7;
    y_exit(i) = elem_B(i)*theta_enter;
    R_o(i) = abs(elem_B(i)*theta_enter);
end

%Transverse magnification
for i = 1:length(z_s)
    if (z_s(i) == f_o)
        trans_mag(i) = (f_tl/f_o)*(z_h/f_1);
    else
        f_e(i) = (f_o*z_s(i))/(f_o-z_s(i));
        f_g(i) = (f_tl*(f_e(i)+d1))/(f_tl-(f_e(i)+d1));
        f_h(i) = (f_1*(f_g(i)+d2))/(f_1-(f_g(i)+d2));
        f_j(i) = (f_d1*(d3+f_h(i)))/(f_d1-(d3+f_h(i)));
        trans_mag(i) = abs((f_e(i)/z_s(i))*(f_g(i)/(f_e(i)+d1))*(f_h(i)/(f_g(i)+d2))*(z_h/(f_h(i)+d3)));
    end
end


%Radius of hologram at camera
y_enter = 0;
theta_enter = NA;
elem_D_term4 = (d1/f_1)*(1-(d2/f_tl));
elem_D_term5 = (d1/f_tl);  
elem_D_term6 = (1-(d2/f_1));
elem_D = zeros(1,length(z_s));
y_exit_camera = zeros(1,length(z_s));
R_spherical_short = zeros(1,length(z_s));
R_spherical_long = zeros(1,length(z_s));
R_h = zeros(1,length(z_s));
OPD_max = zeros(1,length(z_s));
for i = 1: length(z_s) 
    elem_D_term1 = (z_s(i)/f_1)*(1-(d2/f_tl))*(1-(d1/f_o));
    elem_D_term2 = (z_s(i)/f_tl)*(1-(d1/f_o));
    elem_D_term3 = (z_s(i)/f_o)*(1-(d2/f_1));
    elem_D(i) = -elem_D_term1-elem_D_term2-elem_D_term3-elem_D_term4-elem_D_term5+elem_D_term6;
    %Radius of shorter focal length spherical wave at camera
    R_spherical_short(i) = (elem_B(i)*(1-(z_h/f_d1))+(z_h*elem_D(i)))*theta_enter;
    %Radius of shorter focal length spherical wave at camera
    R_spherical_long(i) = (elem_B(i)*(1-(z_h/f_d2))+(z_h*elem_D(i)))*theta_enter;
    R_h(i) = min(abs(R_spherical_short(i)),abs(R_spherical_long(i)));
    OPD_max(i) = sqrt((R_o(i)+R_h(i))^2+z_h^2)- sqrt((R_o(i)-R_h(i))^2+z_h^2);
end

delta_lambda = (690-654)*10^-6;
lambda = 670e-6;
OPD_condition = (lambda^2/delta_lambda);


defocus = (z_s*1e+3-3e+3);
for i = 1: length(z_s)
    magnified_hologram(i) = R_h(i)/trans_mag(i);
    mag_hol(i) = (magnified_hologram(i)*1.0)/sum(magnified_hologram);
end
plot(defocus,trans_mag)

