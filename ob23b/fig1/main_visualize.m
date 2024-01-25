A = 1;
ALPHA = 2*pi;
N = 2;
Y_0 = zeros(N,1);

FIT = @(x) sum(A - A * cos(ALPHA * x) + x.^2, 1);

fig1d = figure;
y = linspace(-1,3,101);
Fy = FIT(y);
plot(y,Fy,'k-');xlabel('$y$');ylabel('$f(y)$');

fig2d = plot_2d(Y_0, FIT, 2, 25);
colormap jet

% fig3d = plot_3d(Y_0, FIT, 2, 70);
% colormap jet

% set(gca,'DataAspectRatio',[1 1 10])
% ExportSetup: custom renderer: painters for high quality PDF
% tightfig.m
% myfigstyle(gcf, 9, 9, 9, 9);