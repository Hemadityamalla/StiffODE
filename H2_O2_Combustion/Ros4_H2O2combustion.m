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
eqNo = 8;
y = zeros(N,eqNo);
facmax = 15;
facmin = 0.5;
safetyFactor = 0.9;
Atol = 1e-6;
Rtol = 1e-3;


%H2-O2 combustion
y(1,:) = [0.3,0.15,0,0,0,0,0,0]; %phi=1 stoichiometric ratio
T = 1e-5;
h = 1e-10;
x = linspace(0,T,N);


%Coefficients for 4th order Rosenbrock Integration(4stage with 3 stage embedded) - GRK4A
%Only for the case of autonomous ODE
gamma = 0.395;
gammaCoeffs = [0, 0, 0, 0;
              -0.7676724, 0, 0, 0;
              -0.8516753, 0.5229673, 0, 0;
              0.28846311, 0.880214e-1, -0.33739, 0];
gammaTilde = gammaCoeffs/gamma;

alpha = [0, 0, 0, 0;
        0.438, 0, 0, 0;
        0.79692045, 0.73076e-1, 0, 0;
        0, 0, 0, 0];
    
c = [0.199293, 0.482645, 0.68061488e-1, 0.25];
chat = [0.3463258, 0.2856931, 0.36798, 0];
order = 4;
%---------------Starting step size prediction-------------%
do = norm(y(1,:));
d1 = 0;
if do < 1e-5 || d1 < 1e-5
    h0 = 1e-6;
else
    h0 = 0.01*(do/d1);
end
startTemp = temperature(x(1));
ytemp = y(1,:) + h0*f(computeRates(startTemp),y(1,:));
midTemp = temperature(x(1)+h0);
d2 = norm(f(computeRates(midTemp),ytemp) - f(computeRates(startTemp),y(1,:)))/h0;
if max(d1,d2) <= 1e-5
    h1 = max(1e-6,h0*1e-3);
else
    p = 4; %Order of the integration scheme being used
    h1 = (0.01/max(d1,d2))^(1/(1+p));
end

h = min(100*h0,h1);


fprintf('Predicted starting step size= %ld \n', h);


%---Start Integration------%
i = 2;
t = 0;
while i < N && t < T
    timeStep = [timeStep,double(h)];
    Sol = zeros(eqNo,1);
    betterSol = zeros(eqNo,1);
    %Computing the solutions for the new step
    Temp = temperature(t); 
    rateConsts = computeRates(Temp);
    J = dfdy(y(i-1,:),rateConsts);
    k = zeros(order,eqNo);
    dx = 1e-5;
    b = zeros(eqNo,1);
    Jt = (f(computeRates(temperature(t+h)),y(i-1,:)) - f(rateConsts,y(i-1,:)))/h;%dfdt(t,y(i-1,:),rateConsts);
    A = eye(length(J)) - h*J*gamma;
    %First Slope computation
    b(:,1) = h*f(rateConsts,y(i-1,:)) + h^2*sum(gammaCoeffs(1,:))*Jt;
    k(1,:) = A\b;
    
    
    %Second Slope computation
    
    b(:,1) = h*f(computeRates(temperature(t+h*alpha(2,1))),y(i-1,:)+alpha(2,1)*k(1,:)) +h^2*sum(gammaCoeffs(2,:))*Jt + gammaTilde(2,1)*k(1,:);
    k(2,:) = A\b - gammaTilde(2,1)*k(1,:)';
    
    
    %Third Slope computation
    b(:,1) = h*f(computeRates(temperature(t + h*(alpha(3,1) + alpha(3,2)))),y(i-1,:)+ alpha(3,1)*k(1,:) +alpha(3,2)*k(2,:)) + h^2*sum(gammaCoeffs(3,:))*Jt + gammaTilde(3,1)*k(1,:)+gammaTilde(3,2)*k(2,:);
    k(3,:) = A\b - gammaTilde(3,1)*k(1,:)'- gammaTilde(3,2)*k(2,:)';
    
    
    %Fourth Slope computation
    b(:,1) = h*f(computeRates(temperature(t + h*(alpha(3,1) + alpha(3,2) + alpha(3,3)))),y(i-1,:) + alpha(3,1)*k(1,:) +alpha(3,2)*k(2,:)) +h^2*sum(gammaCoeffs(4,:))*Jt +  gammaTilde(4,1)*k(1,:)+gammaTilde(4,2)*k(2,:)+gammaTilde(4,3)*k(3,:);
    k(4,:) = A\b - gammaTilde(4,1)*k(1,:)'- gammaTilde(4,2)*k(2,:)' - gammaTilde(4,3)*k(3,:)';
    
    
    %---Actual Integration----
    for ii=1:eqNo
        Sol(ii,1) = y(i-1,ii) + dot(c,k(:,ii));
        betterSol(ii,1) = y(i-1,ii) + dot(chat,k(:,ii));
    end
    %-----Adaptive error analysis-------
    
    sc = Atol*ones(eqNo,1) + max(abs(betterSol(:,1)),abs(y(i-1,:)'))*Rtol;
    
    error = sqrt(sum(((betterSol - Sol)./sc).^2)/eqNo);
    
    if isnan(error)
            fprintf("Solution is Diverging! \n");
            break;
    end
    
    hnew = h*min(facmax,max(facmin,safetyFactor*(1.0/error)^(1/(order+1))));
    if error > 1.0
        fprintf("Rejecting solution and reducing step size! \n");
        hRejected = [hRejected, h];
        tRejected = [tRejected, t];
        facmax = 1;
        h = hnew;
    else
        fprintf("Iteration %d, updating solution! \n",i);
        time = [time, t];
        y(i,:) = betterSol;
        t = t + h;
        i = i + 1;
        hAccepted = [hAccepted, h];
        facmax = 10;
        h = hnew;
    end
    if i >= N
        fprintf("Max iterations exceeded! \n");
        break;
    end
end




fprintf("Final time reached, end of integration!\n");


figure(1)
% plot(time,y(1:length(time),1),'b^-','LineWidth',1,'MarkerSize',3);
% hold on;
% plot(time,y(1:length(time),2),'g^-','LineWidth',1,'MarkerSize',3);
% hold on;
% plot(time,y(1:length(time),3),'r^-','LineWidth',1,'MarkerSize',3);
% hold on;
% plot(time,y(1:length(time),6),'c^-','LineWidth',1,'MarkerSize',3);
% hold on;
plot(time,y(1:length(time),1),'r','LineWidth',2);
hold on;
plot(time,y(1:length(time),2),'g','LineWidth',2);
hold on;
plot(time,y(1:length(time),3),'o','LineWidth',2);
hold on;
plot(time,y(1:length(time),4),'y','LineWidth',2);
hold on;
plot(time,y(1:length(time),5),'c','LineWidth',2);
hold on;
plot(time,y(1:length(time),6),'m','LineWidth',2);
hold on;
plot(time,y(1:length(time),7),'k','LineWidth',2);
hold on;
plot(time,y(1:length(time),8),'b','LineWidth',2);
hold on;
%Solving it with MATLAB solvers
%opts = odeset('RelTol',1e-10,'Abstol',1e-10);
%[t,exact] = ode23s(@f,[0 T],y(1,:),opts);
title('Solution of the Stiff system');
 xlabel('Time t');
 ylabel('Solution y_1');
 %legend('Numerical', 'Exact/MATLAB');
figure(2);
 semilogy(time,hAccepted,'o','MarkerSize',2);
 legend('Accepted');
 hold on;
 semilogy(tRejected,hRejected,'*');
 title('Timesteps- Accepted and rejected')
 xlabel('Iteration');
 ylabel('log(h)');
 legend('Accepted','Rejected');
 
%  figure(3)
%  subplot(4,2,1);
%  loglog(time,y(1:length(time),1),'b^-','LineWidth',1,'MarkerSize',3);
%  hold on;
%  %loglog(t,exact(:,1),'ro','LineWidth',1,'MarkerSize',2);
%  legend('Numerical', 'Exact/MATLAB');
%  subplot(4,2,2);
% loglog(time,y(1:length(time),2),'b^-','LineWidth',1,'MarkerSize',3);
%  hold on;
% %loglog(t,exact(:,2),'ro','LineWidth',1,'MarkerSize',2);
%  legend('Numerical', 'Exact/MATLAB');
%  subplot(4,2,3);
% loglog(time,y(1:length(time),3),'b^-','LineWidth',1,'MarkerSize',3);
%  hold on;
% %loglog(t,exact(:,3),'ro','LineWidth',1,'MarkerSize',2);
%  legend('Numerical', 'Exact/MATLAB');
%  subplot(4,2,4);
% loglog(time,y(1:length(time),4),'b^-','LineWidth',1,'MarkerSize',3);
%  hold on;
% %loglog(t,exact(:,4),'ro','LineWidth',1,'MarkerSize',2);
%  legend('Numerical', 'Exact/MATLAB');
% subplot(4,2,5);
% loglog(time,y(1:length(time),5),'b^-','LineWidth',1,'MarkerSize',3);
%  hold on;
% %loglog(t,exact(:,5),'ro','LineWidth',1,'MarkerSize',2);
%  legend('Numerical', 'Exact/MATLAB');
%  subplot(4,2,6)
% loglog(time,y(1:length(time),6),'b^-','LineWidth',1,'MarkerSize',3);
%  hold on;
% %loglog(t,exact(:,6),'ro','LineWidth',1,'MarkerSize',2);
%  legend('Numerical', 'Exact/MATLAB');
% subplot(4,2,7);
% loglog(time,y(1:length(time),7),'b^-','LineWidth',1,'MarkerSize',3);
%  hold on;
% %loglog(t,exact(:,7),'ro','LineWidth',1,'MarkerSize',2);
%  legend('Numerical', 'Exact/MATLAB');
% subplot(4,2,8);
% loglog(time,y(1:length(time),8),'b^-','LineWidth',1,'MarkerSize',3);
%  hold on;
% %loglog(t,exact(:,8),'ro','LineWidth',1,'MarkerSize',2);
%  legend('Numerical', 'Exact/MATLAB');
