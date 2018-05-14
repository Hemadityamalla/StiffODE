clear;
clc;

N = int64(1e6);
time = [];
%ODE System Init
eqNo = 2;
y = zeros(N,eqNo);
emax = 1e-9;
emin = 1e-10;
hmax = 1e-4;
hmin = 1e-6;
% %Damped- unforced oscillator
% e = 0.1;
% w = 10;
% f = {@(x,y) (y(2));
%      @(x,y) -(2*e*w*y(2) + w^2*y(1))};
% y(1,1) = 0.02;
% y(1,2) = 0;
% T = 4.0;
% h = T/N;
% x = linspace(0,T,N);

%Non-stiff Van-der-Pol Oscillator
mu = 2.0;
f = {@(x,y) (y(2));
     @(x,y) (mu*(1.0 - y(1)^2)*y(2) - y(1))};
y(1,1) = 2.0;
y(1,2) = 0;
T = 20;
h = 1e-10;
x = linspace(0,T,N);

% %Stiff Van-der-Pol Oscillator
% mu = 1000.0;
% f = {@(x,y) (y(2));
%      @(x,y) (mu*(1.0 - y(1)^2)*y(2) - y(1))};
% y(1,1) = 2.0;
% y(1,2) = 0.0;
% z(1,1) = 2.0;
% y(1,2) = 0.0;
% T = 3000;
% h = T/N;
% x = linspace(0,T,N);

%Coefficients for 3rd stage-3rd order Rosenbrock Integration
gamma = 7.88e-1;
alpha = [0.0, 1.0, 1.0];
gammaCoeff = [gamma, -2.11e-1, -1.077e-1];
a = [0,    0, 0;
    1.268, 0, 0;
    1.268, 0, 0];

c = [0,      0,      0;
    -1.6077, 0,      0;
    -3.4641, -1.732, 0];
order = 3;

m = [2.0, 5.7735e-1, 4.2265e-1];
mhat = [2.1132, 1.0, 4.2265e-1];

%---Start Integration------%
i = 2;
t = 0;
while i < N && t < T
    %Solutions used for Rosenbrock33
    Sol = zeros(eqNo,1);
    betterSol = zeros(eqNo,1);
    %--------------
    if h < hmin 
        h = hmin;
    elseif h > hmax
        h = hmax;
    end
    %Computing the solutions for the new step
    J = zeros(eqNo,eqNo);
    for ii = 1:eqNo
        J(ii,:) = jacobianest(f{ii},y(i-1,:),x(i-1));
    end
    k = zeros(order,eqNo);
    dx = 1e-5;
    b = zeros(eqNo,1);
    dfdt = zeros(eqNo,1);
    %First Slope computation
    for ii=1:eqNo
        dfdt(ii,1) = (f{ii}(x(i-1)+dx,y(i-1,:)) - f{ii}(x(i-1),y(i-1,:)))/dx;%Maybe replace this by the derivative function in DERIVSUITE 
        b(ii,1) = f{ii}(x(i-1)+alpha(1)*h,y(i-1,:)) + gammaCoeff(1)*h*dfdt(ii,1);
    end
    k(1,:) = (eye(length(J))/(gamma*h) - J)\b;
    %Second Slope computation
    for ii=1:eqNo
        b(ii,1) = f{ii}(x(i-1)+alpha(2)*h,y(i-1,:)+a(2,1)*k(1,ii)) + gammaCoeff(2)*h*dfdt(ii,1) + (c(2,1)/h)*k(1,ii);
    end
    k(2,:) = (eye(length(J))/(h*gamma) - J)\b;
    %Third Slope computation
    for ii=1:eqNo
        b(ii,1) = f{ii}(x(i-1)+alpha(3)*h,y(i-1,:)+a(2,1)*k(1,ii)) + gammaCoeff(3)*h*dfdt(ii,1) + (c(3,1)*k(1,ii) + c(3,2)*k(2,ii))/h;
    end
    k(3,:) = (eye(length(J))/(h*gamma) - J)\b;
    %---Actual Integration----
    for ii=1:eqNo
        Sol(ii,1) = y(i-1,ii) + dot(m,k(:,ii));
        betterSol(ii,1) = y(i-1,ii) + dot(mhat,k(:,ii));
    end
    error = abs(Sol(1) - betterSol(1));
    if isnan(error)
            fprintf("Solution is Diverging! \n");
            break;
    end
    if error > emax && h > hmin
        h = 0.5*h;
        fprintf("Rejecting solution and reducing step size! \n");
    else
        fprintf("Iteration %d, updating solution! \n",i);
        time = [time, t];
        y(i,:) = betterSol;
        i = i + 1;
        t = t + h;
        if error < emin
            h = 2*h;
            fprintf("Increasing Stepsize! \n");
            
        end
    end
    if i >= N
        fprintf("Max iterations exceeded! \n");
        break;
    end
    %fprintf("Next Iteration. \n");
end




fprintf("Final time reached, end of integration!\n");

%Exact solution
t = linspace(0,T,N);
exact = y(1,1)*exp(-t).*cos(9.95*t);
plot(time,y(1:length(time),1),'-.',t,exact);
hold on;
%[t,exact] = ode23s(@vdp1000,[0 T],[2; 0]);%Check for the value of mu inside @vdp
%plot(t,exact(:,1),'-o');
title('Solution of van der Pol Equation, \mu = 3');
xlabel('Time t');
ylabel('Solution y_1');
legend('Numerical', 'Exact');