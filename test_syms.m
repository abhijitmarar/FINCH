f_o = 3;                            % Focal length of objective (mm)    
NA = 1.42;                          % Numerical aperture of objective
D_bpp = (2*f_o*NA);                 % Diameter of back pupil plane
wave = 670e-6;                      % Wavelength of light (mm)
delta_c = 16e-3;                    % Pixel size of camera (mm)
z_s = 2.990:200e-6:3.010;           % Distance between sample and objective
d_slm = 3;                        % Distance between objective and SLM
f_slm = 300;                        % Focal length of diffractive lens
z_h = 600;                          % Distance between SLM and camera


syms d_SLM f_O z_S z_H f_SLM

f_e = simplify((z_S*f_O)/(f_O-z_S));
num = f_SLM*(f_e+d_SLM);
den = f_SLM-(f_e+d_SLM);
f1 = ((f_SLM*(f_e+d_SLM))/(f_SLM-(f_e+d_SLM)));
f2 = num/den;
num2 = (f1+z_H)*(f_e+d_SLM+z_H);
den2 = f1-f_e-d_SLM;
z_r = num2/den2;


fe_const = double(subs(f_e, [f_O,z_S],[f_o,z_s(1)]));
num_const = subs(num, [f_SLM,f_e, d_SLM],[f_slm,fe_const,d_slm]);
f1_const = subs(f1,[f_SLM,f_e,d_SLM],[f_slm,fe_const,d_slm]);
num2_const = subs(num2, [f1,z_H,f_e,d_SLM],[f1_const,z_h,fe_const,d_slm]);
den2_const = subs(den2, [f1,f_e,d_SLM],[f1_const,fe_const,d_slm]);
z_r_const = double(subs(z_r,[num2,den2],[num2_const,den2_const]));


mag_out_of_focus_num = (f_e*z_H);
mag_out_of_focus_den = (f_e+d_SLM);
mag_out_of_focus = mag_out_of_focus_num/mag_out_of_focus_den;

fe_const = subs(f_e, [f_O,z_S],[f_o,z_s(1)]);
mag_num_cnst = double(subs(mag_out_of_focus_num,[f_e,z_H],[fe_const,z_h]));
mag_den_cnst = z_s(1)*double(subs(mag_out_of_focus_den,[f_e,d_SLM],[fe_const,d_slm]));

