%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CRLB calculation
% 01/14/2020
% Abhijit Marar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIDH with 1 spherical wave, Simple configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEFINE CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_o = 3;                            % Focal length of objective (mm)    
NA = 1.42;                          % Numerical aperture of objective
D_bpp = (2*f_o*NA);                 % Diameter of back pupil plane
wave = 650e-6;                      % Wavelength of light (mm)
k = 2*pi/wave;                      % Wavenumber
z_s = 2.990:50e-6:3.010;           % Distance between sample and objective
d_slm = 3;                          % Distance between objective and SLM
f_slm = 300;                        % Focal length of diffractive lens
z_h = 150;                          % Distance between SLM and camera
r_h = radius_hologram;
alpha = abs(k./(2*z_r));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms r_H z_R K Pi Alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NORMALIZATION CONSTANT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = simplify(1/(r_H^2 + 2*(z_R/K)*sin((K*r_H^2)/(2*z_R))));
normalization_cnst = zeros(1,length(z_s));
for i = 1:length(z_s)
    A_cnst = subs(A,[r_H,z_R,K],[r_h(i),z_r(i),k]);
    normalization_cnst(i) = abs(A_cnst/pi);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sigma_x and Sigma_y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fisher_xx = A*Alpha^2*(((r_H^4)/2)+((1-cos(Alpha*r_H^2))/Alpha^2)-((r_H^2)/Alpha)*sin(Alpha*r_H^2));
F_xx = zeros(1,length(z_s));
sigma_xx = zeros(1,length(z_s));
for i = 1:length(z_s)
    fisher_xx_cnst = subs(fisher_xx, [A,r_H,Alpha],[normalization_cnst(i),r_h(i),alpha(i)]);
    F_xx(i) = 1000*2*pi*abs(fisher_xx_cnst)*trans_mag(i)*trans_mag(i);
    sigma_xx(i) = 1/sqrt(F_xx(i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sigma_z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fisher_zz = (A/z_R^4)*(((r_H^6)/3)-(2*((r_H^2)/Alpha^2)*cos(Alpha*r_H^2))-((((Alpha^2*r_H^4)-2)/(Alpha^3))*sin(Alpha*r_H^2)));
F_zz = zeros(1,length(z_s));
sigma_zz = zeros(1,length(z_s));
for i = 1:length(z_s)
    fisher_zz_cnst = subs(fisher_zz, [A,z_R,r_H,Alpha],[normalization_cnst(i),z_r(i),r_h(i),alpha(i)]);
    F_zz(i) = 1000*((pi*k^2)/2)*abs(fisher_zz_cnst)*df_zr(i)*df_zr(i);
    sigma_zz(i) = 1/sqrt(F_zz(i));
end

defocus = (z_s*1e+3-3e+3);

figure
subplot(1,2,1)
plot(defocus,sigma_xx*1e+6)
axis([-10 10 0 50])
title('CRLB_{xy}')
xlabel('Distance between sample and objective (\mum)')
ylabel('\sigma_{x}, \sigma_{y} (nm)')


subplot(1,2,2)
plot(defocus,sigma_zz*1e+6)
axis([-10 10 0 50])
title('CRLB_{z}')
xlabel('Distance between sample and objective (\mum)')
ylabel('\sigma_{z}(nm)')