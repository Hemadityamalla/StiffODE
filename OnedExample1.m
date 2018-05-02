clear;
clc;


N = 25;
y = zeros(N,1);
y(1) = 0;
h = 1.0/N;
x = linspace(0,1,N);

%Explicit Euler
for i=2:N
    y(i) = y(i-1) + h*f(x(i-1),y(i-1));
end
plot(x,y);
hold on;
