function I = finchH(xs,ys,zs,theta)
% FINCH Microscopy
% Optics Letters 32(8) 912 (2007)

if nargin<=1
    xs = 0;
    ys = 0;
    zs = 0.0;
end
if nargin == 3 || nargin == 0
    theta = 0.0; % radians
end

wave = 0.0005; %wavelength 
a = 430;       %DOE acts as converging diffractive lens with focal length a
f = 250;         %focal length 
d1 = 132;      %distance between positive spherical lens and DOE
d2 = 260;      %distance between SLM and CCD aperture


gamma = (d2-a-zs*(d1*a+d2*f-a*f+d2*a-d1*d2)/f^2)/(1-zs*(a+f-d1)/f^2);
disp(gamma);

N = 2*128;
x = linspace(-1,1,N);
[xx, yy] = meshgrid(x);

amp = exp((1j*pi/wave/gamma)*((xx + a*xs/f).^2 + (yy + a*ys/f).^2)+1j*theta);

I = 2 + amp + conj(amp);   %Eqn(3)

if nargout == 0
    figure();
    imagesc(I);
end