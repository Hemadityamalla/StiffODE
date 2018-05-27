clear;
clc;

x1 = -10;
x2 = -1*x1;

y1 = x1;
y2 = x2;

N = 100;
x = linspace(x1,x2,N);
y = linspace(y1,y2,N);

[x,y] = meshgrid(x,y);
z = x + 1i*y;

stabFun = (1 + (2*z/3) + (z.^2/6))./(1 - z/3);%1 + z + (1/2)*z.^2 + (1/6)*z.^3;

%Plot settings
figure(1)
contour(x,y,abs(stabFun),[0:0.05:1]);
xlim([x1 x2]);
ylim([y1 y2]);
grid on;
xL = xlim;
yL = ylim;
line([0 0], xL,'Color','black');  %x-axis
line(xL, [0 0],'Color','black');  %y-axis