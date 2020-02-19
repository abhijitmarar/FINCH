%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CRLB calculation for SIDH
% 01/14/2020
% Abhijit Marar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEFINE CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_o = 3;                            % Focal length of objective (mm)    
NA = 1.42;                          % Numerical aperture of objective
D_bpp = (2*f_o*NA);                 % Diameter of back pupil plane
wave = 515e-6;                      % Wavelength of light (mm)
k = 2*pi/wave;                      % Wavenumber
z_s = 2.990:200e-6:3.010;           % Distance between sample and objective
d_slm = 3;                          % Distance between objective and SLM
f_slm = 300;                        % Focal length of diffractive lens
z_h = [50,75,100,125,150];          % Distance between SLM and camera before focus
%z_h = [450,475,500,525,550];        % Distance between SLM and camera after focus
r_h = radius_hologram;
r_h(r_h == 0) = eps;
z_r(z_r == 0) = eps;
defocus = (z_s*1e+3-3e+3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms r_H z_R K Pi Alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NORMALIZATION CONSTANT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = simplify(1./((r_H^2)+ 2*(z_R/K)*sin((K*r_H^2)/(2*z_R))));
normalization_cnst = zeros(length(z_h),length(z_s));
alpha = zeros(length(z_h),length(z_s));
for j = 1: length(z_h)
    for i = 1:length(z_s)
        A_cnst = subs(A,[r_H,z_R,K],[r_h(j,i),z_r(j,i),k]);
        normalization_cnst(j,i) = abs(A_cnst/pi);
        alpha(j,i) = abs(k./(2*z_r(j,i)));
    end
end
%normalization_cnst(normalization_cnst == 6.45609149247139e+30) = ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sigma_x and Sigma_y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fisher_xx = A*Alpha^2*(((r_H^4)/2)+((1-cos(Alpha*r_H^2))/Alpha^2)-((r_H^2)/Alpha)*sin(Alpha*r_H^2));
F_xx = zeros(length(z_h),length(z_s));
sigma_xx = zeros(length(z_h),length(z_s));
figure;
p(1) = subplot(1,2,1);
for j = 1:length(z_h)
    for i = 1:length(z_s)
        fisher_xx_cnst = subs(fisher_xx, [A,r_H,Alpha],[normalization_cnst(j,i),r_h(j,i),alpha(j,i)]);
        F_xx(j,i) = 1000*2*pi*abs(fisher_xx_cnst)*trans_mag(j,i)*trans_mag(j,i);
        sigma_xx(j,i) = 1e+6/sqrt(F_xx(j,i));
    end
    plot(defocus,sigma_xx(j,:),'LineWidth',3);
    hold on
end
axis([-10 10 0 50]);
title('CRLB_{xy}');
legend(strcat('z_h =',num2str(z_h'),' mm'),'Location','northwest','FontWeight','bold');
xlabel('Distance between sample and objective (\mum)');
ylabel('\sigma_{x}, \sigma_{y} (nm)');
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sigma_z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fisher_zz = (A/z_R^4)*(((r_H^6)/3)-(2*((r_H^2)/Alpha^2)*cos(Alpha*r_H^2))-((((Alpha^2*r_H^4)-2)/(Alpha^3))*sin(Alpha*r_H^2)));
F_zz = zeros(length(z_h),length(z_s));
sigma_zz = zeros(length(z_h),length(z_s));
p(2) = subplot(1,2,2);
for j = 1: length(z_h)
    for i = 1:length(z_s)
        fisher_zz_cnst = subs(fisher_zz, [A,z_R,r_H,Alpha],[normalization_cnst(j,i),z_r(j,i),r_h(j,i),alpha(j,i)]);
        F_zz(j,i) = 1000*((pi*k^2)/2)*abs(fisher_zz_cnst)*df_zr(j,i)*df_zr(j,i);
        sigma_zz(j,i) = 1e+6/sqrt(F_zz(j,i));
    end
    plot(defocus,sigma_zz(j,:),'LineWidth',3);
    hold on
end
axis([-10 10 0 50])
title('CRLB_{z}')
legend(strcat('z_h =',num2str(z_h'),' mm'),'Location','best','FontWeight','bold');
xlabel('Distance between sample and objective (\mum)')
ylabel('\sigma_{z}(nm)')
