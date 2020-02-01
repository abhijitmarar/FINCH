clear all

%{
SIDH with 2 spherical waves (no tube lens)
Siegel, Nisan, Joseph Rosen, and Gary Brooker. "Reconstruction of objects 
above and below the objective focal plane with dimensional fidelity by FINCH
fluorescence microscopy." Optics express 20.18 (2012): 19822-19835.
%}

f_o = 3;                                % Focal length of objective (mm)    
NA = 1.4;                               % Numerical aperture of objective
R_o = (f_o*NA);
wave = 680e-6;                          % Wavelength of light (mm)
delta_c = 16e-3;                        % Pixel size of camera (mm)
f_d1 = 376.1; f_d2 = 415.7;             % Focal lengths of mirrors
s_fac = (f_d2-f_d1)/(f_d2+f_d1);        % s-factor
z_s = 2.990:200e-6:3.010;
d1 =  50;                                % Distance between objective and interferometer
z_h = (2*f_d1*f_d2)/(f_d1+f_d2);        % Distance between interferometer and camera (maximum overlap)

defocus = (z_s*1e+3-3e+3);

z_r = zeros(1,length(z_s));
trans_mag = zeros(1,length(z_s));
f_e = zeros(1,length(z_s));
f_g = zeros(1,length(z_s));
f_h = zeros(1,length(z_s));
f_j = zeros(1,length(z_s));
gamma = zeros(1,length(z_s));
%Condition to sample finest fringe of Hologram
z_h_min = ((4*R_o*delta_c)/wave);  % Min SLM-camera distance (mm)


%Radius of hologram at Interferometer
y_enter = 0;
theta_enter = NA;
%Element B_o of First final Matrix
y_exit = zeros(1,length(z_s));
for i = 1: length(z_s)
    B_o(i) = z_s(i)*(1-d1/f_o)+d1;
    y_exit(i) = B_o(i)*theta_enter;
    R_o(i) = abs(B_o(i)*theta_enter);
end

plot(defocus,R_o);