f_o = 3;                            % Focal length of objective (mm)    
NA = 1.42;                          % Numerical aperture of objective
D_bpp = (2*f_o*NA);                 % Diameter of back pupil plane
wave = 405e-6;                      % Wavelength of light (mm)
k = 2*pi/wave;                      % Wavenumber

[x,y] = meshgrid(-200:200);
r = sqrt(x.^2 + y.^2);
rm = 500;

%w = rm/10;
term1 = (1 + cos((k*r.^2)/2))*rm;
%term2 = 0.5*tanh((rm - r)/w) + 0.5;
g = term1;
imshow(g,[])