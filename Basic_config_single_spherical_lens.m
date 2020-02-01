%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIDH characterization
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
wave = 515e-6;                      % Wavelength of light (mm)
delta_c = 16e-3;                    % Pixel size of camera (mm)
z_s = 2.990:200e-6:3.010;           % Distance between sample and objective
d_slm = 3;                         % Distance between objective and SLM
f_slm = 300;                        % Focal length of diffractive lens
z_h = 150;                          % Distance between SLM and camera
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RADIUS OF SPHERICAL WAVE AT CAMERA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms d_SLM f_O z_S z_H f_SLM
mat_spherical_wave = simplify([1 z_H;0 1]*[1 0; -1/f_SLM 1]*[1 d_SLM;0 1]*[1 0; -1/f_O 1]*[1 z_S;0 1]);
radius_spherical_wave = zeros(1,length(z_s));
for i = 1:length(z_s)
    mat_spherical_wave_const = subs(mat_spherical_wave,[d_SLM,f_O,z_S,z_H,f_SLM],[d_slm,f_o,z_s(i),z_h,f_slm]);
    radius_spherical_wave(i) = abs(mat_spherical_wave_const(1,2)*NA);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RADIUS OF PLANE WAVE AND HOLOGRAM AT CAMERA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_plane_wave = simplify([1 (d_SLM+z_H);0 1]*[1 0; -1/f_O 1]*[1 z_S;0 1]);
radius_plane_wave = zeros(1,length(z_s));
radius_hologram = zeros(1,length(z_s));
for i = 1:length(z_s)
    mat_plane_wave_const = subs(mat_plane_wave,[d_SLM,f_O,z_S,z_H],[d_slm,f_o,z_s(i),z_h]);
    radius_plane_wave(i) = abs(mat_plane_wave_const(1,2)*NA);
    radius_hologram(i) = min(radius_plane_wave(i),radius_spherical_wave(i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CALCULATION OF RECONSTRUCTION DISTANCE (Z_R)
%Siegel, Nisan, Joseph Rosen, and Gary Brooker. Optics express 20.18 (2012): 19822-19835.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_e = simplify((z_S*f_O)/(f_O-z_S));
f1 = simplify((f_SLM*(f_e+d_SLM))/(f_SLM-(f_e+d_SLM)));
z_r_focus = simplify((f_SLM-z_H));
num2 = (f1+z_H)*(f_e+d_SLM+z_H);
den2 = f1-f_e-d_SLM;
z_r_out_of_focus = num2/den2;
z_r = zeros(1,length(z_s));
for i = 1:length(z_s)
    if (z_s(i) == f_o)
        z_r(i) = double(subs(z_r_focus, [f_SLM,z_H],[f_slm,z_h]));
    else
        fe_const = subs(f_e, [f_O,z_S],[f_o,z_s(i)]);
        f1_const = subs(f1, [f_SLM,f_e,d_SLM],[f_slm,fe_const,d_slm]);
        num2_const = subs(num2, [f1,z_H,f_e,d_SLM],[f1_const,z_h,fe_const,d_slm]);
        den2_const = subs(den2, [f1,f_e,d_SLM],[f1_const,fe_const,d_slm]);
        z_r(i) = double(subs(z_r_out_of_focus, [num2,den2],[num2_const,den2_const]));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CALCULATION OF TRANSVERSE MAGNIFICATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trans_mag = zeros(1,length(z_s));
mag_out_of_focus_num = (f_e*z_H);
mag_out_of_focus_den = (f_e+d_SLM);
for i = 1:length(z_s)
    if (z_s(i) == f_o)
        trans_mag(i) = z_h/f_o;
    else
        fe_const = subs(f_e, [f_O,z_S],[f_o,z_s(i)]);
        mag_num_cnst = double(subs(mag_out_of_focus_num,[f_e,z_H],[fe_const,z_h]));
        mag_den_cnst = z_s(i)*double(subs(mag_out_of_focus_den,[f_e,d_SLM],[fe_const,d_slm]));
        trans_mag(i) = mag_num_cnst/mag_den_cnst;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defocus = (z_s*1e+3-3e+3);

h(1) = subplot(3,2,1);
plot(defocus,radius_spherical_wave)
title('Radius of spherical wave')
xlabel('Distance between sample and objective (\mum)')
ylabel('Radius(mm)')


h(2) = subplot(3,2,2);
plot(defocus,radius_plane_wave)
title('Radius of plane wave')
xlabel('Distance between sample and objective (\mum)')
ylabel('Radius(mm)')

h(3) = subplot(3,2,3);
plot(defocus,radius_hologram)
title('Radius of hologram')
xlabel('Distance between sample and objective (\mum)')
ylabel('Radius(mm)')

h(4) = subplot(3,2,4);
plot(defocus,trans_mag)
title('Transverse magnification')
xlabel('Distance between sample and objective (\mum)')
ylabel('M_{T}')

h(5) = subplot(3,2,5);
plot(defocus,z_r)
title('Reconstruction Distance')
xlabel('Distance between sample and objective (\mum)')
ylabel('z_{r}(mm)')
pos = get(h,'Position');
new = mean(cellfun(@(v)v(1),pos(1:2)));
set(h(5),'Position',[new,pos{end}(2:end)])