%{
function I = finchH_objective(x_s,y_s,z_s,theta)
% FINCH Microscopy
% Optics Express 19(6) 5047 (2011)

if nargin<=1
    x_s = 0;
    y_s = 0;
    z_s = 0;
end
if nargin == 3 || nargin == 0
    theta = 0.0; % radians
end
%}

wave = 0.0005; %wavelength 
f_d = 400;     %focal length of difractive spherical lens (fd = 2*zh is optimal)
f_o = 3.000;   %focal length of objective
d1 = 300;      %distance between objective and SLM
z_h = 200;     %distance between SLM and CCD aperture
x_s = 0.0;
y_s = 0.0;
z_s = 3.000;
theta = 0.0;
z_r = 0;       %reconstruction distance
trans_mag = 0; %transverse magnification

%{
if (z_s > f_o)
    f_e = (z_s*f_o)/(f_o-z_s);
elseif (z_s < f_0)
    f_e = (z_s*f_o)/abs(f_o+z_s);
end
%}
%z_s =2.995;

if (abs(z_s) == 3.000)
    %Simplification f_o = z_s
    z_r = (z_h - f_d);
    disp(z_r);
    trans_mag = z_h/f_o;
else
    %Equation 9
    f_e = (z_s*f_o)/(f_o-z_s);
    f_1 = (f_d*(f_e+d1))/(f_d-(f_e+d1));
    z_r = -(((f_1+z_h)*(f_e+d1+z_h))/(f_1-f_e-d1));
    disp(z_r)
    trans_mag = (z_h*f_e)/(z_s*(f_e+d1));
end

N = 2*256;
x = linspace(-1,1,N);
[xx, yy] = meshgrid(x);

amp = exp((1j*pi/wave/z_r)*((xx - trans_mag*x_s).^2 + (yy - trans_mag*y_s).^2)+1j*theta);

I = 2 + amp + conj(amp);   %Eqn(8) simplified
if nargout == 0
    figure();
    imagesc(I);

end
