%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CRLB calculation for a Gaussian PSF with noise
% 02/07/2020
% Abhijit Marar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

f_o = 3;                                % Focal length of objective (mm)    
NA = 1.42;                              % Numerical aperture of objective
wave = 515e-6;                          % Wavelength of light (mm)
FWHM = (0.61*wave)/NA;                  % FWHM
sigma_o = FWHM/(2*sqrt(2*log(2)));      % standard deviation of Gaussian
z = 2.999:20e-6:3.001;                  % Distance between sample and objective
z_s = (z-3.000);
N = 1000;                               % No. of photons                                
d = 800e-6;                             % Depth of focus (mm)
bg = 200;                               % No. of bg photons/area    


%{
sigma_g = sigma_o;
limit = 6*sigma_g;
q_gauss = @(x,y) (1/(2*pi*sigma_g^2)).*exp(-(x.^2+y.^2)/(2.*sigma_g.^2));
dq_dx =  @(x,y) (-x./(2.*pi.*sigma_g.^4)).*exp(-(x.^2+y.^2)/(2.*sigma_g.^2)); 
F_xx = @(x,y)((N./(bg+q_gauss(x,y))).*dq_dx(x,y).*dq_dx(x,y));
check = integral2(F_xx,-limit,limit,-limit,limit);
%check = integral2(F_xx,-0.1,0.1,-0.1,0.1);
crlb = 1e+6/sqrt(check); % CRLB in nm
%}





sigma_g = zeros(1,length(z_s));
fisher_x = zeros(1,length(z_s));
CRLB_x = zeros(1,length(z_s));
for i = 1:length(z_s)
    sigma_g(i) = sigma_o*sqrt(1+(z_s(i)/d)^2);
    q_gauss = @(x,y) (1./(2*pi*sigma_g(i).^2)).*exp(-(x.^2+y.^2)/(2*sigma_g(i).^2));
    dq_dx = @(x,y) (-x./(2.*pi.*sigma_g(i).^4)).*exp(-(x.^2+y.^2)/(2.*sigma_g(i).^2));
    F_xx =  @(x,y) ((N./(bg+q_gauss(x,y))).*dq_dx(x,y).*dq_dx(x,y));
    fisher_x(i) = integral2(F_xx,-6*sigma_g(i),6*sigma_g(i),-6*sigma_g(i),6*sigma_g(i));
    CRLB_x(i) = 1/sqrt(fisher_x(i));
end



%{
defocus = z_s*1.e+3;

figure
subplot(1,2,1)
plot(defocus,CRLB_x.*1e+6)
%axis([-1 1 0 30])
title('CRLB_{xy}')
xlabel('Distance between sample and objective (\mum)')
ylabel('\sigma_{x}, \sigma_{y} (nm)')

subplot(1,2,2)
plot(defocus,sigma_z*1e+6)
axis([-1 1 0 100])
title('CRLB_{z}')
xlabel('Distance between sample and objective (\mum)')
ylabel('\sigma_{z}(nm)')
%}

