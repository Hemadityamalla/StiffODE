clear;
clc;
h = 0.01;
tfinal = 1000.0;
N = int64(tfinal/h);
y = zeros(N+1,2);
t = linspace(0,tfinal,N+1);
%Initial conditions
y(1,1) = 2;
y(1,2) = 0;

%Initialize system
[y(2,1),y(2,2)] = EulerIntegration([y(1,1),y(1,2)],h,t(1));

%Main integration
for i=1:N-1
    initialGuess = [y(i+1,1),y(i+1,2)];
    parameters = [h, initialGuess, t(i+2),1000]';
    [y(i+2,1),y(i+2,2)] = newtonRaphson(initialGuess, parameters);
end

%Plotting solution
plot(t,y(:,1),'.')%,t,y(:,2));
% hold on;
% time = linspace(0,tfinal,N);
% exact = y(1,1)*exp(-time).*cos(9.95*time);
% plot(time,exact);

%Necessary functions
function [fy1,fy2] = f(t,y)

    %fy1 = y(1) - y(2);
    %fy2 = 4*y(1) - 3*y(2);
    %fy1 = -2*y(2)^3;
    %fy2 = 2*y(1) - y(2)^3;
    %Damped Oscillator
%     e = 0.1;
%     w = 10;
%     fy1 = y(2);
%     fy2 = -(2*e*w*y(2) + w^2*y(1));
    %van der Pol
     mu = 1000.0;
     fy1 = y(2);
     fy2 = (mu*(1.0 - y(1)^2)*y(2) - y(1));


end

function [y1,y2] = EulerIntegration(yn, h,t)
    [fy1,fy2] = f(t,yn);
    %y1 = yn(1) + h*fy1;
    %y2 = yn(2) + h*fy2;
    %Second order
    [fy1,fy2] = f(t+0.5*h,yn+0.5*h*[fy1,fy2]);
    y1 = yn(1) + h*fy1;
    y2 = yn(2) + h*fy2;
end

function [x,y] = newtonRaphson(guess, params)
    h = params(1,1);
    xn = params(2,1);
    yn = params(3,1);
    t = params(4,1);
    mu = params(5,1);
    TOL = 1e-8;
    root = [0,0];
    iter = 0;
    [fx,fy] = f(t,[guess(1),guess(2)]);
    F = -[xn,yn] - h*[fx,fy] + [guess(1),guess(2)]; %This will be zero at extrema- so take care of that
    while abs(max(F)) > TOL || abs(max(root-guess)) > TOL
        if iter ~=0
            guess = root;
        end
        %J = [1,6*h*guess(2)^2;-2*h,1+3*h*guess(2)^2];
        %J = [1,-h;h*100,1+2*h]; %Jacobian for the damped oscillator
        J = eye(length(guess))- h*[0, 1; -2*mu*guess(2)*guess(1), mu*(1 - guess(1)^2)];
        root = [guess(1),guess(2)] - (J\F')';
        iter = iter+1;
        [fx,fy]= f(t,[guess(1),guess(2)]);
        F = -[xn,yn] - h*[fx,fy] + [root(1),root(2)];
    end
    %fprintf('Number of iterations: %d\n',iter);
    x = root(1);
    y = root(2);
end