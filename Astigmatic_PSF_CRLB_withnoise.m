%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CRLB calculation for Astigmatic PSF with noise
% 02/07/2020
% Abhijit Marar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

f_o = 3;                                % Focal length of objective (mm)    
NA = 1.42;                              % Numerical aperture of objective
wave = 670e-6;                          % Wavelength of light (mm)
FWHM = (0.61*wave)/NA;                  % FWHM
sigma_o = FWHM/(2*sqrt(2*log(2)));      % standard deviation of Gaussian
z = 2.990:50e-6:3.010;                  % Distance between sample and objective
z_s = (z-3.000);
N = 1000;                               % No. of photons                                
d = 800e-6;                             % Depth of focus (mm)
gamma = 400e-6;                         % Amount of astigmatism (mm)
bg = 200;                               % No. of bg photons/area

sigma_x = zeros(1,length(z_s));
sigma_y = zeros(1,length(z_s));
dsigmax_dz = zeros(1,length(z_s));
dsigmay_dz = zeros(1,length(z_s));
CRLB_x = zeros(1,length(z_s));
CRLB_y = zeros(1,length(z_s));
CRLB_z = zeros(1,length(z_s));
epsilon = zeros(1,length(z_s));
fisher_x = zeros(1,length(z_s));
fisher_y = zeros(1,length(z_s));
fisher_z = zeros(1,length(z_s));

for i = 1:length(z_s)
    sigma_x(i) = sigma_o*sqrt(1+((z_s(i)+gamma)/d)^2);
    sigma_y(i) = sigma_o*sqrt(1+((z_s(i)-gamma)/d)^2);
    epsilon(i) = sqrt(sigma_y(i)/sigma_x(i));
    q_astigmatic = @(x,y) (1./(2*pi*sigma_x(i).*sigma_y(i))).*exp(-(((x.^2)/(2*sigma_x(i).^2))+((y.^2)/(2*sigma_y(i).^2))));
    dq_dx = @(x,y) (-x./(2.*pi.*sigma_y(i).*sigma_x(i).^3)).*exp(-(((x.^2)/(2*sigma_x(i).^2))+((y.^2)/(2*sigma_y(i).^2))));
    dq_dy = @(x,y) (-y./(2.*pi.*sigma_x(i).*sigma_y(i).^3)).*exp(-(((x.^2)/(2*sigma_x(i).^2))+((y.^2)/(2*sigma_y(i).^2))));
    dq_dsigmax = @(x,y) q_astigmatic(x,y).*((x.^2/sigma_x(i).^3)-(1/sigma_x(i)));
    dq_dsigmay = @(x,y) q_astigmatic(x,y).*((y.^2/sigma_y(i).^3)-(1/sigma_y(i)));
    dsigmax_dz(i) = ((sigma_o*(z_s(i)+gamma))/d^2)*(1+((z_s(i)+gamma)/d)^2)^(-1/2);
    dsigmay_dz(i) = ((sigma_o*(z_s(i)-gamma))/d^2)*(1+((z_s(i)-gamma)/d)^2)^(-1/2);
    dq_dz = @(x,y) (dq_dsigmax(x,y).*dsigmax_dz(i)) + (dq_dsigmay(x,y).*dsigmay_dz(i));
    F_xx =  @(x,y) ((N./(bg+q_astigmatic(x,y))).*dq_dx(x,y).*dq_dx(x,y));
    F_yy =  @(x,y) ((N./(bg+q_astigmatic(x,y))).*dq_dy(x,y).*dq_dy(x,y));
    F_zz =  @(x,y) ((N./(bg+q_astigmatic(x,y))).*dq_dz(x,y).*dq_dz(x,y));
    fisher_x(i) = integral2(F_xx,-6*sigma_x(i),6*sigma_x(i),-6*sigma_y(i),6*sigma_y(i));
    fisher_y(i) = integral2(F_yy,-6*sigma_x(i),6*sigma_x(i),-6*sigma_y(i),6*sigma_y(i));
    fisher_z(i) = integral2(F_zz,-6*sigma_x(i),6*sigma_x(i),-6*sigma_y(i),6*sigma_y(i));
    CRLB_x(i) = 1/sqrt(fisher_x(i));
    CRLB_y(i) = 1/sqrt(fisher_y(i));
    CRLB_z(i) = 1/sqrt(fisher_z(i));
end

defocus = z_s*1.e+3;

figure
subplot(1,2,1)
plot(defocus,CRLB_x.*1e+6,'LineWidth',3)
hold on
plot(defocus,CRLB_y.*1e+6,'LineWidth',3)

axis([-2 2 0 30])
title('CRLB_{xy}')
xlabel('Distance between sample and objective (\mum)')
ylabel('\sigma_{x}, \sigma_{y} (nm)')
hold off

subplot(1,2,2)
plot(defocus,CRLB_z*1e+6,'LineWidth',3)
axis([-2 2 0 100])
title('CRLB_{z}')
xlabel('Distance between sample and objective (\mum)')
ylabel('\sigma_{z}(nm)')
    
    
    
    
    
    
    