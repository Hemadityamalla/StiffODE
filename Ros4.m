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
safetyFactor = 0.9;
Atol = 1e-8;
Rtol = 1e-6;
Dtol = 1e-4;

% %Damped- unforced oscillator
% e = 0.1;
% w = 10;
% f = {@(x,y) (y(2));
%      @(x,y) -(2*e*w*y(2) + w^2*y(1))};
% y(1,1) = 0.02;
% y(1,2) = 0;
% T = 4.0;
% h = 1e-4;
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
T = 2000;
h = 1e-3;
x = linspace(0,T,N);


%Coefficients for 4th order Rosenbrock Integration(4stage with 3 stage embedded) - GRK4A
%Only for the case of autonomous ODE
gamma = 0.395;
gammaCoeff = [0, 0, 0, 0;
              -0.7676724, 0, 0, 0;
              -0.8516753, 0.5229673, 0, 0;
              0.28846311, 0.880214e-1, -0.33739, 0];
gammaTilde = gammaCoeff/gamma;

alpha = [0, 0, 0, 0;
        0.438, 0, 0, 0;
        0.79692045, 0.73076e-1, 0, 0;
        0, 0, 0, 0];
    
c = [0.199293, 0.482645, 0.68061488e-1, 0.25];
chat = [0.3463258, 0.2856931, 0.36798, 0];
order = 4;

%---Start Integration------%
i = 2;
t = 0;
while i < N && t < T
    timeStep = [timeStep,double(h)];
    Sol = zeros(eqNo,1);
    betterSol = zeros(eqNo,1);
    %Computing the solutions for the new step
%     J = zeros(eqNo,eqNo);
%     for ii = 1:eqNo
%         J(ii,:) = jacobianest(f{ii},y(i-1,:),x(i-1));
%     end
    J = [0, 1; -2*mu*y(i-1,2)*y(i-1,1), mu*(1 - y(i-1,1)^2)];
%    J = [0,1;-w^2, -2*e*w];
    k = zeros(order,eqNo);
    dx = 1e-5;
    b = zeros(eqNo,1);
    dfdt = zeros(eqNo,1);
    A = eye(length(J)) - h*J*gamma;
    %First Slope computation
    for ii=1:eqNo
        dfdt(ii,1) = (f{ii}(x(i-1)+dx,y(i-1,:)) - f{ii}(x(i-1),y(i-1,:)))/dx;%Maybe replace this by the derivative function in DERIVSUITE 
        b(ii,1) = h*f{ii}(x(i-1),y(i-1,:));
    end
    %c1 = [c1, cond(A)];
    k(1,:) = A\b;
    %Second Slope computation
    for ii=1:eqNo
        b(ii,1) = h*f{ii}(x(i-1),y(i-1,:)+alpha(2,1)*k(1,ii)) + gammaTilde(2,1)*k(1,ii);
    end
    k(2,:) = A\b - gammaTilde(2,1)*k(1,:)';
    %Third Slope computation
    for ii=1:eqNo
        b(ii,1) = h*f{ii}(x(i-1),y(i-1,:) + alpha(3,1)*k(1,ii) +alpha(3,2)*k(2,ii))+gammaTilde(3,1)*k(1,ii)+gammaTilde(3,2)*k(2,ii);
    end
    k(3,:) = A\b - gammaTilde(3,1)*k(1,:)'- gammaTilde(3,2)*k(2,:)';
    %Fourth Slope computation
    for ii=1:eqNo
        b(ii,1) = h*f{ii}(x(i-1),y(i-1,:) + alpha(3,1)*k(1,ii) +alpha(3,2)*k(2,ii))+gammaTilde(4,1)*k(1,ii)+gammaTilde(4,2)*k(2,ii)+gammaTilde(4,3)*k(3,ii);
    end
    k(4,:) = A\b - gammaTilde(4,1)*k(1,:)'- gammaTilde(4,2)*k(2,:)' - gammaTilde(4,3)*k(3,:)';
    %---Actual Integration----
    for ii=1:eqNo
        Sol(ii,1) = y(i-1,ii) + dot(c,k(:,ii));
        betterSol(ii,1) = y(i-1,ii) + dot(chat,k(:,ii));
    end
    %-----Adaptive error analysis-------
    
    sc = [Atol + max(abs(betterSol(1,1)),abs(y(i-1,1)))*Rtol+max(abs(f{1}(x(i),betterSol(:,1))),abs(f{1}(x(i-1),y(i-1,:))))*Dtol,...
          Atol + max(abs(betterSol(2,1)),abs(y(i-1,2)))*Rtol+max(abs(f{2}(x(i),betterSol(:,1))),abs(f{2}(x(i-1),y(i-1,:))))*Dtol];

%     sc = [Atol + max(abs(betterSol(1,1)),abs(y(i-1,1)))*Rtol,...
%           Atol + max(abs(betterSol(2,1)),abs(y(i-1,2)))*Rtol];
    
    error = sqrt(sum(((betterSol - Sol)./sc').^2)/eqNo);
    
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

%Exact solution
t = linspace(0,T,N);
exact = y(1,1)*exp(-t).*cos(9.95*t);
figure(2)
subplot(2,1,1);
plot(time,y(1:length(time),1));%,'-o',t,exact);
hold on;
[t,exact] = ode23s(@vdp1000,[0 T],[2; 0]);%Check for the value of mu inside @vdp
plot(t,exact(:,1),'-o')%,t,exact(:,2),'^',t,exact(:,3),'*');
title('Solution of van der Pol Equation, \mu = 1000');
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

% %Plotting the condition numbers
% figure(3)
% grid on;
% plot(c1,'ro-');
% hold on;
% plot(timeStep*1e2,'b*-')
% % hold on;
% % loglog(vpa(c2),'b*');
% % hold on;
% % loglog(vpa(c3),'g^');
% xlabel('Iteration');
% ylabel('2-norm condition number');
% legend('Condition number','scaled Stepsize');



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
