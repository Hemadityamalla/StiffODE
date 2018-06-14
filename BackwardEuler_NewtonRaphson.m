clear;
clc;

h = 0.1;
tfinal = 5.0;
n = int64(tfinal/h);

y = zeros(n,1);
x=  zeros(n,1);
t = linspace(0,tfinal,n);

%Initial Conditions
y(1,1) = 0;
x(1,1) = 0.02;

for i=1:n-1
    %Solving for the i+1th values using newton-raphson
    initialGuess = [x(i,1),y(i,1)];
    parameters = [h, initialGuess, t(i+1)]';
    [x(i+1,1),y(i+1,1)] = newtonRaphson(initialGuess, parameters);
end

plot(t,x,t,y);

function [fx,fy] = f(t,x,y)
    %fx = -2*y^3;
    %fy = 2*x - y^3;
    e = 0.1;
    w = 10;
    fx = y;
    fy = -(2*e*w*y + w^2*x);
end

function [x,y] = newtonRaphson(guess, params)
    h = params(1,1);
    xn = params(2,1);
    yn = params(3,1);
    t = params(4,1);
    TOL = 1e-8;
    root = [0,0];
    iter = 0;
    [fx,fy] = f(t,guess(1),guess(2));
    F = -[xn,yn] - h*[fx,fy] + [guess(1),guess(2)];
    while abs(max(F)) > TOL
        if iter ~=0
            guess = root;
        end
       % J = [1,6*h*guess(2)^2;-2*h,1+3*h*guess(2)^2];
       J = [0,1;-100,-2];
        root = [guess(1),guess(2)] - (J\F')';
        iter = iter+1;
        [fx,fy]= f(t,guess(1),guess(2));
        F = -[xn,yn] - h*[fx,fy] + [root(1),root(2)];
    end
    x = root(1);
    y = root(2);
end

