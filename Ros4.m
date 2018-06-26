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
Rtol = 1e-4;
Dtol = 1e-4;

% %Damped- unforced oscillator
% e = 0.1;
% w = 10;
% y(1,1) = 0.02;
% y(1,2) = 0;
% T = 4.0;
% h = 1e-4;
% x = linspace(0,T,N);

% %Non-stiff Van-der-Pol Oscillator
% mu = 2.0;
% y(1,1) = 2.0;
% y(1,2) = 0;
% T = 20;
% h = 1e-5;
% x = linspace(0,T,N);

% 
% % Stiff Van-der-Pol Oscillator
% mu = 1000;
% y(1,1) = 2.0;
% y(1,2) = 0;
% T = 3000;
% h = 1e-3;
% x = linspace(0,T,N);

%Stiff HIRES
y(1,:) = [1,0,0,0,0,0,0,0.0057];
T = 321.8122;
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
%---------------Starting step size prediction-------------%
do = norm(y(1,:));
d1 = 0;
if do < 1e-5 || d1 < 1e-5
    h0 = 1e-6;
else
    h0 = 0.01*(do/d1);
end
ytemp = y(1,:) + h0*fHIRES(x(1),y(1,:));
d2 = norm(fHIRES(x(1)+h0,ytemp) - fHIRES(x(1),y(1,:)))/h0;
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

    %J = [0, 1; -2*mu*y(i-1,2)*y(i-1,1), mu*(1 - y(i-1,1)^2)];
    %J = [0,1;-w^2, -2*e*w];
    J = [-1.71,  0.43, 8.32, 0, 0, 0, 0, 0;...
          1.71, -8.75,    0, 0, 0, 0, 0, 0;...
             0,     0,-10.03, 0.43,0.035,0,0,0;...
             0,8.32,1.71,-1.12,0,0,0,0;...
             0,0,0,0,-1.745,0.43,0.43,0;...
             0,0,0,0.69,1.71,-0.43-280*y(i-1,8),0.69,-280*y(i-1,6);...
             0,0,0,0,0,280*y(i-1,8),-1.81,280*y(i-1,6);...
             0,0,0,0,0,-280*y(i-1,8),1.81,-280*y(i-1,6)];
    k = zeros(order,eqNo);
    dx = 1e-5;
    b = zeros(eqNo,1);
    dfdt = zeros(eqNo,1);
    A = eye(length(J)) - h*J*gamma;
    %First Slope computation
    %b(:,1) = h*fDO(x(i-1),y(i-1,:),e,w);
    %b(:,1) = h*fVDP(x(i-1),y(i-1,:),mu);
    b(:,1) = h*fHIRES(x(i-1),y(i-1,:));
    k(1,:) = A\b;
    
    
    %Second Slope computation
    %b(:,1) = h*fDO(x(i-1),y(i-1,:)+alpha(2,1)*k(1,:),e,w) + gammaTilde(2,1)*k(1,:);
    %b(:,1) = h*fVDP(x(i-1),y(i-1,:)+alpha(2,1)*k(1,:),mu) + gammaTilde(2,1)*k(1,:);
    b(:,1) = h*fHIRES(x(i-1),y(i-1,:)+alpha(2,1)*k(1,:)) + gammaTilde(2,1)*k(1,:);
    k(2,:) = A\b - gammaTilde(2,1)*k(1,:)';
    
    
    %Third Slope computation
    %b(:,1) = h*fDO(x(i-1),y(i-1,:)+ alpha(3,1)*k(1,:) +alpha(3,2)*k(2,:),e,w)+gammaTilde(3,1)*k(1,:)+gammaTilde(3,2)*k(2,:);
    %b(:,1) = h*fVDP(x(i-1),y(i-1,:)+ alpha(3,1)*k(1,:) +alpha(3,2)*k(2,:),mu)+gammaTilde(3,1)*k(1,:)+gammaTilde(3,2)*k(2,:);
    b(:,1) = h*fHIRES(x(i-1),y(i-1,:)+ alpha(3,1)*k(1,:) +alpha(3,2)*k(2,:))+gammaTilde(3,1)*k(1,:)+gammaTilde(3,2)*k(2,:);
    k(3,:) = A\b - gammaTilde(3,1)*k(1,:)'- gammaTilde(3,2)*k(2,:)';
    
    
    %Fourth Slope computation
    %b(:,1) = h*fDO(x(i-1),y(i-1,:) + alpha(3,1)*k(1,:) +alpha(3,2)*k(2,:),e,w)+gammaTilde(4,1)*k(1,:)+gammaTilde(4,2)*k(2,:)+gammaTilde(4,3)*k(3,:);
    %b(:,1) = h*fVDP(x(i-1),y(i-1,:) + alpha(3,1)*k(1,:) +alpha(3,2)*k(2,:),mu)+gammaTilde(4,1)*k(1,:)+gammaTilde(4,2)*k(2,:)+gammaTilde(4,3)*k(3,:);
    b(:,1) = h*fHIRES(x(i-1),y(i-1,:) + alpha(3,1)*k(1,:) +alpha(3,2)*k(2,:))+gammaTilde(4,1)*k(1,:)+gammaTilde(4,2)*k(2,:)+gammaTilde(4,3)*k(3,:);
    k(4,:) = A\b - gammaTilde(4,1)*k(1,:)'- gammaTilde(4,2)*k(2,:)' - gammaTilde(4,3)*k(3,:)';
    
    
    %---Actual Integration----
    for ii=1:eqNo
        Sol(ii,1) = y(i-1,ii) + dot(c,k(:,ii));
        betterSol(ii,1) = y(i-1,ii) + dot(chat,k(:,ii));
    end
    %-----Adaptive error analysis-------
    
%     sc = [Atol + max(abs(betterSol(1,1)),abs(y(i-1,1)))*Rtol+max(abs(f{1}(x(i),betterSol(:,1))),abs(f{1}(x(i-1),y(i-1,:))))*Dtol,...
%           Atol + max(abs(betterSol(2,1)),abs(y(i-1,2)))*Rtol+max(abs(f{2}(x(i),betterSol(:,1))),abs(f{2}(x(i-1),y(i-1,:))))*Dtol];

    %sc = zeros(eqNo,1);
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
        %fprintf("Iteration %d, updating solution! \n",i);
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
t = linspace(0,T,1000);
exact = y(1,1)*exp(-t).*cos(9.95*t);
figure(2)
subplot(2,1,1);
plot(time,y(1:length(time),1),'b^-','LineWidth',1,'MarkerSize',3);
hold on;
% plot(t,exact,'r.');
% hold on;
opts = odeset('RelTol',1e-10,'Abstol',1e-10);
[t,exact] = ode23s(@fHIRESMATLAB,[0 T],y(1,:),opts);
plot(t,exact(:,1),'ro','LineWidth',1,'MarkerSize',2);
title('Solution of the HIRES system');
 xlabel('Time t');
 ylabel('Solution y_1');
 legend('Numerical', 'Exact/MATLAB');
subplot(2,1,2);
 semilogy(time,hAccepted,'o','MarkerSize',2);
 legend('Accepted');
 hold on;
 semilogy(tRejected,hRejected,'*');
 title('Timesteps- Accepted and rejected')
 xlabel('Iteration');
 ylabel('log(h)');
 %semilogy(t);
 legend('Accepted','Rejected');
 
%  figure(3)
%  subplot(4,2,1);
%  loglog(time,y(1:length(time),1),'b^-','LineWidth',1,'MarkerSize',3);
%  hold on;
%  loglog(t,exact(:,1),'ro','LineWidth',1,'MarkerSize',2);
%  legend('Numerical', 'Exact/MATLAB');
%  subplot(4,2,2);
% loglog(time,y(1:length(time),2),'b^-','LineWidth',1,'MarkerSize',3);
%  hold on;
% loglog(t,exact(:,2),'ro','LineWidth',1,'MarkerSize',2);
%  legend('Numerical', 'Exact/MATLAB');
%  subplot(4,2,3);
% loglog(time,y(1:length(time),3),'b^-','LineWidth',1,'MarkerSize',3);
%  hold on;
% loglog(t,exact(:,3),'ro','LineWidth',1,'MarkerSize',2);
%  legend('Numerical', 'Exact/MATLAB');
%  subplot(4,2,4);
% loglog(time,y(1:length(time),4),'b^-','LineWidth',1,'MarkerSize',3);
%  hold on;
% loglog(t,exact(:,4),'ro','LineWidth',1,'MarkerSize',2);
%  legend('Numerical', 'Exact/MATLAB');
% subplot(4,2,5);
% loglog(time,y(1:length(time),5),'b^-','LineWidth',1,'MarkerSize',3);
%  hold on;
% loglog(t,exact(:,5),'ro','LineWidth',1,'MarkerSize',2);
%  legend('Numerical', 'Exact/MATLAB');
%  subplot(4,2,6)
% loglog(time,y(1:length(time),6),'b^-','LineWidth',1,'MarkerSize',3);
%  hold on;
% loglog(t,exact(:,6),'ro','LineWidth',1,'MarkerSize',2);
%  legend('Numerical', 'Exact/MATLAB');
% subplot(4,2,7);
% loglog(time,y(1:length(time),7),'b^-','LineWidth',1,'MarkerSize',3);
%  hold on;
% loglog(t,exact(:,7),'ro','LineWidth',1,'MarkerSize',2);
%  legend('Numerical', 'Exact/MATLAB');
% subplot(4,2,8);
% loglog(time,y(1:length(time),8),'b^-','LineWidth',1,'MarkerSize',3);
%  hold on;
% loglog(t,exact(:,8),'ro','LineWidth',1,'MarkerSize',2);
%  legend('Numerical', 'Exact/MATLAB');



 
 
 
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

function yDO = fDO(x,y,e,w)

yDO = [y(2), -(2*e*w*y(2) + w^2*y(1))];

end

function yVDP = fVDP(x,y,mu)

yVDP = [y(2), (mu*(1.0 - y(1)^2)*y(2) - y(1))];

end

function yHIRES = fHIRES(x,y)

yHIRES = [-1.71*y(1)+0.43*y(2)+8.32*y(3)+0.0007,...
           1.71*y(1)-8.75*y(2),...
          -10.03*y(3)+0.43*y(4)+0.035*y(5),...
          8.32*y(2)+1.71*y(3)-1.12*y(4),...
          -1.745*y(5)+0.43*y(6)+0.43*y(7),...
          -280*y(6)*y(8)+0.69*y(4)+1.71*y(5)-0.43*y(6)+0.69*y(7),...
          280*y(6)*y(8)-1.81*y(7),...
          -280*y(6)*y(8)+1.81*y(7)];
 

end

function yHIRES = fHIRESMATLAB(x,y)

yHIRES = [-1.71*y(1)+0.43*y(2)+8.32*y(3)+0.0007,...
           1.71*y(1)-8.75*y(2),...
          -10.03*y(3)+0.43*y(4)+0.035*y(5),...
          8.32*y(2)+1.71*y(3)-1.12*y(4),...
          -1.745*y(5)+0.43*y(6)+0.43*y(7),...
          -280*y(6)*y(8)+0.69*y(4)+1.71*y(5)-0.43*y(6)+0.69*y(7),...
          280*y(6)*y(8)-1.81*y(7),...
          -280*y(6)*y(8)+1.81*y(7)]';
 

end