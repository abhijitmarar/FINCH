%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CRLB calculation for SIDH with Poisson noise
% 02/12/2020
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
z_h = 150;                          % Distance between SLM and camera
r_h = radius_hologram;
N = 10000;                          % No. of photons in hologram
bg = 5;                             % Bg photons            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = zeros(1,length(z_s));
fisher_x = zeros(1,length(z_s));
CRLB_x = zeros(1,length(z_s));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(z_s)
    A(i) = 1/(pi*(r_h(i)^2+2*(z_r(i)/k)*sin((k/(2*z_r(i)))*r_h(i)^2)));
    q_sidh = @(a,b) A(i).*(1+cos((k/(2*z_r(i))).*(a.^2+b.^2)));
    dq_dxs = @(a,b) 2.*A(i).*trans_mag(i).*(k/(2*z_r(i))).*sin((k/(2*z_r(i))).*(a.^2+b.^2)).*a;
    F_xx = @(a,b) ((N./(q_sidh(a,b)+bg)).*dq_dxs(a,b).*dq_dxs(a,b));
    polar_Fxx = @(theta,r) F_xx(r.*cos(theta),r.*sin(theta)).*r;
    fisher_x(i) = integral2(polar_Fxx,0,2*pi,0,r_h(i));
    CRLB_x(i) = 1.e+6/sqrt(fisher_x(i));
end

defocus = (z_s*1e+3-3e+3);

figure
subplot(1,2,1)
plot(defocus,CRLB_x)
axis([-10 10 10 50])
title('CRLB_{xy}')
xlabel('Distance between sample and objective (\mum)')
ylabel('\sigma_{x}, \sigma_{y} (nm)')

%{
subplot(1,2,2)
plot(defocus,CRLB_z*1e+6)
axis([-1 1 0 100])
title('CRLB_{z}')
xlabel('Distance between sample and objective (\mum)')
ylabel('\sigma_{z}(nm)')
%}
