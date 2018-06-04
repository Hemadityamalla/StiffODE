clear;
clc;

x1 = -20;
x2 = -1*x1;

y1 = x1;
y2 = x2;

N = 100;
x = linspace(x1,x2,N);
y = linspace(y1,y2,N);

[x,y] = meshgrid(x,y);
z = x + 1i*y;

gamma = 0.25;%1 + 1/sqrt(2.0);

stabFun = (1 + (1 - 2*gamma)*z + (0.5 - 2*gamma + gamma^2)*z.^2)./(1 - gamma*z).^2;%1 + z + (1/2)*z.^2 + (1/6)*z.^3;

%Plot settings
figure(1)
contour(x,y,abs(stabFun),[0:0.01:1]);
xlim([x1 x2]);
ylim([y1 y2]);
grid on;
xL = xlim;
yL = ylim;
line([0 0], xL,'Color','black');  %x-axis
line(xL, [0 0],'Color','black');  %y-axis