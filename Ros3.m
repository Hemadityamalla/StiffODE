clear;
clc;

N = int64(1e6);
time = [];
c1 = [];
c2 = [];
c3 = [];
timeStep = [];
%ODE System Init
eqNo = 2;
y = zeros(N,eqNo);
emax = 1e-6;
emin = 1e-7;
hmax = 0.2;
hmin = 1e-7;
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

% %Non-stiff Van-der-Pol Oscillator
% mu = 2.0;
% f = {@(x,y) (y(2));
%      @(x,y) (mu*(1.0 - y(1)^2)*y(2) - y(1))};
% y(1,1) = 2.0;
% y(1,2) = 0;
% T = 20;
% h = 1e-5;
% x = linspace(0,T,N);
% 
%Scaled - Stiff Van-der-Pol Oscillator
mu = 1e6;
f = {@(x,y) (y(2));
     @(x,y) mu*((1.0 - y(1)^2)*y(2) - y(1))};
y(1,1) = 2.0;
y(1,2) = -0.66;
T = 2.5;
h = T/N;
x = linspace(0,T,N);

%Stiff Chemical Kinetics Reaction
% f = {@(x,y) (-0.04*y(1) + 1e4*y(2)*(3));
%     @(x,y) (0.04*y(1) - 1e4*y(2)*y(3) - 3e7*y(2)^2);
%     @(x,y) (3e7*y(2)^2)
%     };
% eqNo = 3;
% y = zeros(N,eqNo);
% y(1,1) = 1;
% y(1,2) = 0;
% y(1,3) = 0;
% T = 10;
% h = 1e-2;
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
    timeStep = [timeStep,double(h)];
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
%     J = zeros(eqNo,eqNo);
%     for ii = 1:eqNo
%         J(ii,:) = jacobianest(f{ii},y(i-1,:),x(i-1));
%     end
    J = [0, 1; -2*mu*y(i-1,2)*y(i-1,1), mu*(1 - y(i-1,1)^2)];
    k = zeros(order,eqNo);
    dx = 1e-5;
    b = zeros(eqNo,1);
    dfdt = zeros(eqNo,1);
    %First Slope computation
    for ii=1:eqNo
        dfdt(ii,1) = (f{ii}(x(i-1)+dx,y(i-1,:)) - f{ii}(x(i-1),y(i-1,:)))/dx;%Maybe replace this by the derivative function in DERIVSUITE 
        b(ii,1) = f{ii}(x(i-1)+alpha(1)*h,y(i-1,:)) + gammaCoeff(1)*h*dfdt(ii,1);
    end
    A = (eye(length(J))/(gamma*h) - J);
    c1 = [c1, cond(A)];
    k(1,:) = A\b;
    %Second Slope computation
    for ii=1:eqNo
        b(ii,1) = f{ii}(x(i-1)+alpha(2)*h,y(i-1,:)+a(2,1)*k(1,ii)) + gammaCoeff(2)*h*dfdt(ii,1) + (c(2,1)/h)*k(1,ii);
    end
    A = (eye(length(J))/(gamma*h) - J);
    %c2 = [c2, vpa(cond(A))];
    k(2,:) = A\b;
    %Third Slope computation
    for ii=1:eqNo
        b(ii,1) = f{ii}(x(i-1)+alpha(3)*h,y(i-1,:)+a(2,1)*k(1,ii)) + gammaCoeff(3)*h*dfdt(ii,1) + (c(3,1)*k(1,ii) + c(3,2)*k(2,ii))/h;
    end
    A = (eye(length(J))/(gamma*h) - J);
    %c3 = [c3, vpa(cond(A))];
    k(3,:) = A\b;
    %---Actual Integration----
    for ii=1:eqNo
        Sol(ii,1) = y(i-1,ii) + dot(m,k(:,ii));
        betterSol(ii,1) = y(i-1,ii) + dot(mhat,k(:,ii));
    end
    %-----Adaptive error analysis-------
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
%t = linspace(0,T,N);
%exact = y(1,1)*exp(-t).*cos(9.95*t);
figure(1)
plot(time,y(1:length(time),1));%,'-.',t,exact);
hold on;
[t,exact] = ode23s(@vdp1000,[0 T],[2; -0.66]);%Check for the value of mu inside @vdp
plot(t,exact(:,1),'-o')%,t,exact(:,2),'^',t,exact(:,3),'*');
title('Solution of van der Pol Equation, \mu = 1000');
xlabel('Time t');
ylabel('Solution y_1');
legend('Numerical', 'Exact/MATLAB');
% 
%Plotting the condition numbers
figure(2)
grid on;
plot(c1,'ro-');
hold on;
plot(timeStep*1e2,'b*-')
% hold on;
% loglog(vpa(c2),'b*');
% hold on;
% loglog(vpa(c3),'g^');
xlabel('Iteration');
ylabel('2-norm condition number');
legend('Condition number','scaled Stepsize');



% %Writing condition numbers to files
% fileID = fopen('C1_DO_1e6.txt','w');
% fprintf(fileID,'%6s %12s\n','iteration','condNo_k1');
% fprintf(fileID,'%6.2f %12.8f\n',[linspace(0,length(c1),length(c1));vpa(c1)]);
% fclose(fileID);
% 
% fileID = fopen('C2_DO_1e6.txt','w');
% fprintf(fileID,'%6s %12s\n','iteration','condNo_k2');
% fprintf(fileID,'%6.2f %12.8f\n',[linspace(0,length(c2),length(c2));vpa(c2)]);
% fclose(fileID);
% 
% fileID = fopen('C3_DO_1e6.txt','w');
% fprintf(fileID,'%6s %12s\n','iteration','condNo_k3');
% fprintf(fileID,'%6.2f %12.8f\n',[linspace(0,length(c3),length(c3));vpa(c3)]);
% fclose(fileID);

