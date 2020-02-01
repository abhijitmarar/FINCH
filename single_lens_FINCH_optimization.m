clear all

% SIDH with 1 spherical wave

f_o = 3;                            % Focal length of objective (mm)    
NA = 1.42;                           % Numerical aperture of objective
D_bpp = (2*f_o*NA);                 % Diameter of back pupil plane
wave = 670e-6;                      % Wavelength of light (mm)
delta_c = 16e-3;                    % Pixel size of camera (mm)
f_1 = 300;                          % Focal lengths of mirrors
f_tl = 180;                         % Focal length of tube lens
f_4 = 120;                           % Focal length of L4
z_s = 2.990:200e-6:3.010;
d1 = 183;                           % Distance between objective and tube lens
d2 = 300;                           % Distance between tube lens and L4
d3 = 120;                           % Distance between L4 and Interferometer
z_h = 150;                          % Distance between interferometer and camera (maximum overlap)
%z_h = 200;
z_r = zeros(1,length(z_s));

f_e = zeros(1,length(z_s));
f_g = zeros(1,length(z_s));
f_h = zeros(1,length(z_s));
f_k = zeros(1,length(z_s));


%{
%Condition to sample finest fringe of Hologram
z_h_min = ((4*D_bpp*(*delta_c)/wave);  % Min SLM-camera distance (mm)
%}

%Reconstruction Distance
for i = 1:length(z_s)
    if (z_s(i) == f_o)
        z_r(i) = -(z_h-f_1);
    else
        f_e(i) = (f_o*z_s(i))/(f_o-z_s(i));
        f_g(i) = (f_tl*(f_e(i)+d1))/(f_tl-(f_e(i)+d1));
        f_h(i) = (f_4*(f_g(i)+d2))/(f_4-(f_g(i)+d2));
        f_k(i) = (f_1*(f_h(i)+d3))/(f_1-(f_h(i)+d3));
        z_r(i) = ((z_h+f_k(i))*(f_h(i)+d3+z_h))/(f_k(i)-f_h(i)-d3);
    end
end

%Radius of hologram at Interferometer
y_enter = 0;
theta_enter = zeros(1,length(z_s));
%Element B of Jones Matrix
elem_B_term5 = (d1*(1-(d3/f_4))*(1-(d2/f_tl)));
elem_B_term6 = (d3*(1-(d1/f_tl)));
elem_B_term7 = (d2*(1-(d3/f_4)));
elem_B = zeros(1,length(z_s));
y_exit = zeros(1,length(z_s));
R_o = zeros(1,length(z_s));
for i = 1: length(z_s)
    elem_B_term1 = (1-(d1/f_o))*(1-(d3/f_4))*(1-(d2/f_tl))*(z_s(i)); 
    elem_B_term2 = (d3/f_tl)*(1-(d1/f_o))*(z_s(i));
    elem_B_term3 = (d2/f_o)*(1-(d3/f_4))*(z_s(i));
    elem_B_term4 = (d3/f_o)*(z_s(i));  
    elem_B(i) = elem_B_term1-elem_B_term2-elem_B_term3-elem_B_term4+elem_B_term5+elem_B_term6+elem_B_term7;
    theta_enter(i) = (f_o*NA)/z_s(i);
    y_exit(i) = elem_B(i)*theta_enter(i);
    R_o(i) = abs(elem_B(i)*theta_enter(i));
end



%Radius of hologram at camera
y_enter = 0;
%theta_enter = NA;
theta_enter = zeros(1,length(z_s));
elem_D_term4 = (d1/f_4)*(1-(d2/f_tl));
elem_D_term5 = (d1/f_tl);  
elem_D_term6 = (1-(d2/f_4));
elem_D = zeros(1,length(z_s));
y_exit_camera = zeros(1,length(z_s));
R_spherical = zeros(1,length(z_s));
R_h = zeros(1,length(z_s));
OPD_max = zeros(1,length(z_s));
for i = 1: length(z_s) 
    elem_D_term1 = (z_s(i)/f_4)*(1-(d2/f_tl))*(1-(d1/f_o));
    elem_D_term2 = (z_s(i)/f_tl)*(1-(d1/f_o));
    elem_D_term3 = (z_s(i)/f_o)*(1-(d2/f_4));
    elem_D(i) = -elem_D_term1-elem_D_term2-elem_D_term3-elem_D_term4-elem_D_term5+elem_D_term6;
    theta_enter(i) = (f_o*NA)/z_s(i);
    %Radius of plane wave at camera
    y_exit_camera(i) = (elem_B(i)+(z_h*elem_D(i)))*theta_enter(i);
    %Radius of spherical wave at camera
    R_spherical(i) = (elem_B(i)*(1-(z_h/f_1))+(z_h*elem_D(i)))*theta_enter(i);
    %R_h(i) = abs(min(R_spherical(i),y_exit_camera(i)));
    R_h(i) = min(abs(R_spherical(i)),abs(y_exit_camera(i)));
    OPD_max(i) = sqrt((R_o(i)+R_h(i))^2+z_h^2)- sqrt((R_o(i)-R_h(i))^2+z_h^2);
end
delta_lambda = (690-654)*10^-6;
lambda = 670e-6;
OPD_condition = (lambda^2/delta_lambda);

%Tranverse Magnification
trans_mag = zeros(1,length(z_s));
for i = 1:length(z_s)
    if (z_s(i) == f_o)
        trans_mag(i) = (f_tl/f_o)*(z_h/f_4);
    else
        f_e(i) = (f_o*z_s(i))/(f_o-z_s(i));
        f_g(i) = (f_tl*(f_e(i)+d1))/(f_tl-(f_e(i)+d1));
        f_h(i) = (f_4*(f_g(i)+d2))/(f_4-(f_g(i)+d2));
        trans_mag(i) = abs((z_h*f_e(i)*f_g(i)*f_h(i))/(z_s(i)*(d1+f_e(i))*(d2+f_g(i))*(d3+f_h(i))));
    end
end

defocus = (z_s*1e+3-3e+3);
figure
plot(defocus,uint8(trans_mag))
figure
plot(defocus,R_h)
%figure
%plot(defocus,trans_mag)
%xlabel('Defocus (um)')
%ylabel('Transverse Magnification(mm)')





