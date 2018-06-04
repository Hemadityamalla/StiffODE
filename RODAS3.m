clear;
clc;

N = 40000;%int64(4e4);
time = [];
%ODE System Init
eqNo = 2;
y = zeros(N,eqNo);
% %Damped- unforced oscillator
% e = 0.1;
% w = 10;
% f = {@(x,y) (y(2));
%      @(x,y) -(2*e*w*y(2) + w^2*y(1))};
% y(1,1) = 0.02;
% y(1,2) = 0;
% T = 4.0;
% h = 1e-8;
% x = linspace(0,T,N);

%Stiff Van-der-Pol Oscillator
mu = 1000.0;
f = {@(x,y) (y(2));
     @(x,y) (mu*(1.0 - y(1)^2)*y(2) - y(1))};
y(1,1) = 2.0;
y(1,2) = 0.0;
z(1,1) = 2.0;
y(1,2) = 0.0;
T = 3000;
h = 1e-10;
x = linspace(0,T,N);

gamma = 0.5;
b = [5/6,-1/6,-1/6,1/2];


%Start integration
t = 0;
for i = 2:N
    %Solutions used for Rosenbrock12
    secondOrderSol = zeros(eqNo,1);

    %---Using the Rosenbrock2 method----
    order = 2;
    J = [0, 1; -2*mu*y(i-1,2)*y(i-1,1), mu*(1 - y(i-1,1)^2)];
    %J = [0,1;-w^2, -2*e*w];
        %Computing the slopes ki
    f1 = [f{1}(x(i-1),y(i-1,:)), f{2}(x(i-1), y(i-1,:))]';
    k1 = (eye(length(J)) - gamma*h*J)\f1;
    k2 = (eye(length(J)) - gamma*h*J)\(f1 + J*k1);
    f2 = [f{1}(x(i-1)+h,y(i-1,:)+h*k1(1)),...
         f{2}(x(i-1)+h, y(i-1,:)+h*k1(2))]';
    k3 = (eye(length(J)) - gamma*h*J)\(f2 - 0.25*J*k1 - 0.25*J*k2);
    f3 = [f{1}(x(i-1)+h,y(i-1,:)+h*(0.75*k1-0.25*k2+0.5*k3)),...
          f{2}(x(i-1)+h, y(i-1,:)+h*(0.75*k1-0.25*k2+0.5*k3))]';
    k4 = (eye(length(J)) - gamma*h*J)\(f3 + (J/12)*(k1+k2-8*k3));
    %---Actual Integration----
    for ii=1:eqNo
        secondOrderSol(ii,1) = y(i-1,ii) + h*dot(b,[k1(ii),k2(ii),k3(ii),k4(ii)]);
    end
    if max(abs(secondOrderSol)) > 2
        fprintf("Solution is diverging\n");
        break;
    end
    y(i,:) = secondOrderSol(:,1);
    if i >= N
        fprintf("Max iterations exceeded! \n");
        break;
    end
    fprintf("Iteration %d. \n",i);
    t = t + h;
end
fprintf("Final time reached, end of integration!\n");

%Exact solution
t = linspace(0,T,N);
exact = y(1,1)*exp(-t).*cos(9.95*t);
plot(t,y(:,1))%,'-.',t,exact);
hold on;
% [t,y] = ode15s(@vdp,[0 T],[2; 0]);%Check for the value of mu inside @vdp
% plot(t,y(:,1),'-o');
% title('Solution of van der Pol Equation, \mu = 1000');
xlabel('Time t');
ylabel('Solution y_1');
legend('Numerical', 'Exact');