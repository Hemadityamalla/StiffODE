N = 40000;%int64(4e4);
time = [];
%ODE System Init
eqNo = 2;
y = zeros(N,eqNo);
%Damped- unforced oscillator
e = 0.1;
w = 10;
f = {@(x,y) (y(2));
     @(x,y) -(2*e*w*y(2) + w^2*y(1))};
y(1,1) = 0.02;
y(1,2) = 0;
T = 4.0;
h = T/N;
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


%Start integration
t = 0;
for i = 2:N
    %Solutions used for Rosenbrock12
    secondOrderSol = zeros(eqNo,1);

    %---Using the Rosenbrock3 method----
    order = 3;
    
        %Computing the Jacobian
    J = zeros(eqNo,eqNo);
    for ii = 1:eqNo
        J(ii,:) = jacobianest(f{ii},y(i-1,:),x(i-1));
    end
        %Computing the slopes ki
    kr = zeros(order,eqNo);
    gamma = 7.886751345948129e-01; %chosen such that the scheme is stable for larger step sizes.
    dx = 1e-5; %Used to estimate the derivative
    b = zeros(eqNo,1);
    dfdt = zeros(eqNo,1);
        %First slope- k1
        gamma1 = gamma;
    for ii=1:eqNo
        dfdt(ii,1) = (f{ii}(x(i-1)+dx,y(i-1,:)) - f{ii}(x(i-1),y(i-1,:)))/dx;%Maybe replace this by the derivative function in DERIVSUITE 
        b(ii,1) = f{ii}(x(i-1),y(i-1,:)) + h*gamma1*dfdt(ii,1);
    end
    kr(1,:) = (eye(length(J))/(h*gamma) - J)/b';
        %Second slope- k2
        alpha21 = 1.267949192431123e0;
        c21 = -1.607695154586736e+0;
        gamma2 = -2.113248654051871e-01;
    for ii=1:eqNo
        b(ii,1) = f{ii}(x(i-1)+h,y(i-1,:)+alpha21*kr(1,ii)) + (c21/h)*kr(1,ii) + h*gamma2*dfdt(ii,1); %Extra term arises for the non-autonomous case.
    end
    kr(2,:) = (eye(length(J))/(h*gamma) - J)/b';
        %Third slope- k3
        alpha31 = 1.267949192431123e0;
        alpha32 = 0.000000000000000e+0;
        c31 = -1.607695154586736e+0;
        c32 = -1.732050807568877e+0;
        gamma3 = -1.077350269189626e+0;
    for ii=1:eqNo
        b(ii,1) = f{ii}(x(i-1)+h,y(i-1,:)+alpha31*kr(1,ii)+alpha32*kr(2,ii)) + (c31/h)*kr(1,ii) + (c32/h)*kr(2,ii) + h*gamma3*dfdt(ii,1); %Extra term arises for the non-autonomous case.
    end
    kr(3,:) = (eye(length(J))/(h*gamma) - J)/b';
    
    %---Actual Integration----
    m1 = 2.0;
    m2 = 5.773502691896258e-01;
    m3 = 4.226497308103742e-01;
    for ii=1:eqNo
        secondOrderSol(ii,1) = y(i-1,ii) + m1*kr(1,ii)+m2*kr(2,ii)+m3*kr(3,ii);
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
plot(t,y(:,1),'-.',t,exact);
% hold on;
% [t,y] = ode15s(@vdp1000,[0 T],[2; 0]);
% plot(t,y(:,1),'-o');
% title('Solution of van der Pol Equation, \mu = 1000');
xlabel('Time t');
ylabel('Solution y_1');
legend('Numerical', 'Exact');