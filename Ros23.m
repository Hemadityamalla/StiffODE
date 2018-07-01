clear;
clc;


N = int64(1e6);
time = [];
c1 = [];
c2 = [];
c3 = [];
hAccepted = [];
hRejected = [];
timeStep = [];
tRejected = [];
%ODE System Init
eqNo = 2;
y = zeros(N,eqNo);
facmax = 15;
facmin = 0.5;
safetyFactor = 0.8;
Atol = 1e-6;
Rtol = 1e-3;

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

% %Non-stiff Van-der-Pol Oscillator
% mu = 2.0;
% f = {@(x,y) (y(2));
%      @(x,y) (mu*(1.0 - y(1)^2)*y(2) - y(1))};
% y(1,1) = 2.0;
% y(1,2) = 0;
% T = 20;
% h = 1e-5;
% x = linspace(0,T,N);

% Stiff Van-der-Pol Oscillator
mu = 1000;
f = {@(x,y) (y(2));
     @(x,y) mu*(1-y(1)^2)*y(2) - y(1)};
y(1,1) = 2.0;
y(1,2) = 0;
T = 3000;
h = 1e-3;
x = linspace(0,T,N);


d = 1/(2 + sqrt(2));
e32 = 6 + sqrt(2);
order = 3;

%---Start Integration------%
i = 2;
t = 0;
while i < N && t < T
    timeStep = [timeStep,double(h)];
    Sol = zeros(eqNo,1);
    betterSol = zeros(eqNo,1);
    %Computing the solutions for the new step

    J = [0, 1; -2*mu*y(i-1,2)*y(i-1,1), mu*(1 - y(i-1,1)^2)];
    %J = [0,1;-w^2, -2*e*w];
    k = zeros(order,eqNo);
    F = zeros(order,eqNo);
    dx = 1e-5;
    %b = zeros(eqNo,1);
    dfdt = zeros(eqNo,1);
    W = eye(length(J)) - h*J*d;
    %First Slope computation
    for ii=1:eqNo
        dfdt(ii,1) = (f{ii}(x(i-1)+dx,y(i-1,:)) - f{ii}(x(i-1),y(i-1,:)))/dx;%Maybe replace this by the derivative function in DERIVSUITE 
        F(1,ii) = f{ii}(x(i-1),y(i-1,:));
    end
    %c1 = [c1, cond(A)];
    k(1,:) = W\F(1,:)';
    %Second Slope computation
    for ii=1:eqNo
        F(2,ii) = f{ii}(x(i-1),y(i-1,:)+0.5*h*k(1,:));
    end
    k(2,:) = W\(F(2,:)' - k(1,:)') + k(1,:)';
    %---Actual Integration----
    for ii=1:eqNo
        Sol(ii,1) = y(i-1,ii) + h*k(2,ii);
    end
    %Third Slope computation
    for ii=1:eqNo
        F(3,ii) = f{ii}(x(i),Sol(:,1));
    end
    k3 = W\(F(3,:)' - e32*(k(2,:) - F(2,:)') - 2*(k(1,:) - F(1,:)'));
    %-----Adaptive error analysis-------
    
    threshold = Atol/Rtol;
    error = (h/6)*norm(k(1,:) - 2*k(2,:) + k(3,:)) / max(max(norm(Sol),norm(y(i-1,:))),threshold);
    
    if isnan(error)
            fprintf("Solution is Diverging! \n");
            break;
    end
    
    %hnew = h*min(facmax,max(facmin,safetyFactor*(Rtol/error)^(1/(2+1))));
    hnew = max(h*facmin,h*max(0.5,0.8*(Rtol/error)^(1/(2+1))));
    if error > Rtol
        fprintf("Rejecting solution and reducing step size! \n");
        hRejected = [hRejected, h];
        tRejected = [tRejected, t];
        facmax = 1;
        h = hnew;
    else %Updating the accepted solution
        fprintf("Iteration %d, updating solution! \n",i);
        time = [time, t];
        y(i,:) = Sol;
        t = t + h;
        i = i + 1;
        hAccepted = [hAccepted, h];
        h = hnew;
    end
    if i >= N
        fprintf("Max iterations exceeded! \n");
        break;
    end
end




fprintf("Final time reached, end of integration!\n");

%Exact solution
% t = linspace(0,T,N);
% exact = y(1,1)*exp(-t).*cos(9.95*t);
figure(2)
subplot(2,1,1);
plot(time,y(1:length(time),1),'b^-','LineWidth',1,'MarkerSize',3);
hold on;
opts = odeset('RelTol',1e-12,'Abstol',1e-14);
%[tym,exact] = ode23s(@vdpstiff,[0 T],[2; 0],opts);%Check for the value of mu inside @vdp
A = csvread("VDP_Exact.dat");
plot(A(:,1),A(:,2),'r-')%,t,exact(:,2),'^',t,exact(:,3),'*');
title('Solution of Stiff Equation, \mu = 1000');
 xlabel('Time t');
 ylabel('Solution y_1');
 legend('Numerical', 'Exact/MATLAB');
subplot(2,1,2);
 semilogy(time,hAccepted,'o');
 legend('Accepted');
 hold on;
 semilogy(tRejected,hRejected,'*');
 title('Timesteps- Accepted and rejected')
 xlabel('Iteration');
 ylabel('log(h)');
 %semilogy(t);
 legend('Accepted','Rejected');

function dydt = vdpns(t,y)
% Evaluate the van der Pol ODEs for mu = 2.
 
dydt = [y(2); 2*(1-y(1)^2)*y(2) - y(1)];
end

function dydt = vdpstiff(t,y)
% Evaluate the van der Pol ODEs for mu = 1000.
 
dydt = [y(2); 1000*(1-y(1)^2)*y(2) - y(1)];
end