function I = finchH(xs,ys,zs,theta)
% FINCH Microscopy
% Optics Letters 32(8) 912 (2007)

if nargin<=1
    xs = 0;
    ys = 0;
    zs = 0.0;
end
if nargin == 3 || nargin == 0
    theta = 0.5; % radians
end

wave = 0.0005;
a = 200; 
f = 3;


d1 = f;
d2 = 180;

gamma = (d2-a-zs*(d1*a+d2*f-a*f+d2*a-d1*d2)/f^2)/(1-zs*(a+f-d1)/f^2);

N = 128;
x = linspace(-0.2,0.2,N);
[xx, yy] = meshgrid(x);

amp = exp((1j*pi/wave/gamma)*((xx + a*xs/f).^2 + (yy + a*ys/f).^2)+1j*theta);

I = 2 + amp + conj(amp);

if nargout == 0
    figure();
    imagesc(I);
end