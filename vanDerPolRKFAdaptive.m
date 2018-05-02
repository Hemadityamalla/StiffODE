clear;
clc;


N = int64(1e6);
time = [];
%ODE System Init
eqNo = 2;
y = zeros(N,eqNo);
emax = 1e-7;
emin = 1e-17;
hmax = 1.0;
hmin = 1e-17;
% %Damped- unforced oscillator
% e = 0.1;
% w = 10;
% f = {@(x,y) (y(2));
%      @(x,y) -(2*e*w*y(2) + w^2*y(1))};
% y(1,1) = 0.02;
% y(1,2) = 0;
% T = 4.0;
% h = 1e-2;
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

% %Stiff Chemical Kinetics Reaction
% f = {@(x,y) (-0.04*y(1) + 1e4*y(2)*(3));
%     @(x,y) (0.04*y(1) - 1e4*y(2)*y(3) - 3e7*y(2)^2);
%     @(x,y) (3e7*y(2)^2)
%     };
% eqNo = 3;
% y = zeros(N,eqNo);
% y(1,1) = 1;
% y(1,2) = 0;
% y(1,3) = 0;
% T = 1e2;
% h = 1e-2;
% x = linspace(0,T,N);
 
 



%RKFehlberg Data
order = 6; %This is not necessarily the order of accuracy, it is just the number of slopes needed
c = [0.0,1.0/4.0,3.0/8.0,12.0/13.0,1.0,0.5];
b1 = [25.0/216.0,1408.0/2565.0,2197.0/4101.0,-1.0/5.0];
b2 = [16.0/135.0,6656.0/12825.0,28561.0/56430.0,-9.0/50.0,2.0/55.0];
k = zeros(order,eqNo);
A = [0,0,0,0,0,0;
    0.25,0,0,0,0,0;
    3.0/32.0,9.0/32.0,0,0,0,0;
    1932.0/2197.0,-7200.0/2197.0,7296.0/2197.0,0,0,0;
    439.0/216.0,-8.0,3680.0/513.0,-845.0/4104.0,0,0;
    -8.0/27.0,2.0,-3544.0/2565.0,1859.0/4104.0,-11.0/40.0,0];


%Start integration
i = 2;
t = 0;
while i < N && t < T
    fourthOrderSol = zeros(eqNo,1);
    fifthOrderSol = zeros(eqNo,1);
    %Solutions used for Rosenbrock12
    secondOrderSol = zeros(eqNo,1);
    firstOrderSol = zeros(eqNo,1);
    %--------------
    if h < hmin 
        h = hmin;
    elseif h > hmax
        h = hmax;
    end
    %Computing the solutions for the new step
    %----Using the RKF45 method---
%     for ii=1:eqNo
%             k(1,ii) = f{ii}(x(i-1),y(i-1,:));
%     end
%     for j = 2:order
%         temp = repmat(A(j,1:j-1),eqNo,1);
%         for ii=1:eqNo
%             k(j,ii) = f{ii}(x(i-1) + c(j)*h, y(i-1,:) + h*dot(temp',k(1:j-1,:),1));
%         end
%     end
    %---Using the Rosenbrock12 method----
        %Computing the Jacobian
    J = zeros(eqNo,eqNo);
    for ii = 1:eqNo
        J(ii,:) = jacobianest(f{ii},y(i-1,:),x(texi-1));
    end
    kr = zeros(2,eqNo);
    gamma = 1.0 + 1/sqrt(2.0);
    dx = 1e-5;
    b = zeros(eqNo,1);
    dfdt = zeros(eqNo,1);
    for ii=1:eqNo
        dfdt(ii,1) = (f{ii}(x(i-1)+dx,y(i-1,:)) - f{ii}(x(i-1),y(i-1,:)))/dx;%Maybe replace this by the derivative function in DERIVSUITE 
        b(ii,1) = f{ii}(x(i-1),y(i-1,:)) + gamma*h*dfdt(ii,1);
    end
    kr(1,:) = J/b';
    for ii=1:eqNo
        b(ii,1) = f{ii}(x(i-1)+h,y(i-1,:)+h*kr(1,ii)) - gamma*h*dfdt(ii,1) - 2*kr(1,ii);
    end
    kr(2,:) = J/b';
    
    %---Actual Integration----
    for ii=1:eqNo
        fourthOrderSol(ii,1) = y(i-1,ii) + h*(dot(b1,[k(1,ii),k(3,ii),k(4,ii),k(5,ii)]));
        fifthOrderSol(ii,1) = y(i-1,ii) + h*(dot(b2,[k(1,ii),k(3,ii),k(4,ii),k(5,ii),k(6,ii)]));
        secondOrderSol(ii,1) = y(i-1,ii) + (0.5*h)*(3*kr(1,ii) + k(2,ii));
        firstOrderSol(ii,1) = y(i-1,ii) + h*kr(1,ii);
    end
    error = abs(secondOrderSol(1) - firstOrderSol(1));
    %error = abs(fourthOrderSol(1) - fifthOrderSol(1));
    if isnan(error)
            fprintf("Solution is Diverging! \n");
            break;
    end
    if error > emax && h > hmin
        h = 0.5*h;
        fprintf("Rejecting solution and reducing step size! \n");
    else
        fprintf("Iteration %d, updating solution! \n",i);
        %y(i,:) = fifthOrderSol;
        time = [time, t];
        y(i,:) = secondOrderSol;
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


plot(time,y(1:length(time),1),'-.');