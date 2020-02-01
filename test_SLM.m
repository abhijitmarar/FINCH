x = linspace(-3,3,N);
[xx, yy] = meshgrid(x);
n_rows = 512;
n_col = 512;
init = zeros(n_rows,n_col);

wave = 0.515;
a = 400;
theta1 = 0;


R_phase1 = ((pi/wave/a)*(xx.^2+yy.^2)-theta1);
imagesc(R_phase1)