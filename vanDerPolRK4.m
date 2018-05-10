clear;
clc;


N = 4000;
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

% %Non-stiff Van-der-Pol Oscillator
% mu = 1.0;
% f = {@(x,y) (y(2));
%      @(x,y) (mu*(1.0 - y(1)^2)*y(2) - y(1))};
% y(1,1) = 2.0;
% y(1,2) = 0;
% T = 20;
% h = T/N;
% x = linspace(0,T,N);

% %Stiff Van-der-Pol Oscillator
% mu = 1000.0;
% f = {@(x,y) (y(2));
%      @(x,y) (mu*(1.0 - y(1)^2)*y(2) - y(1))};
% y(1,1) = 2.0;
% y(1,2) = 0;
% T = 3000;
% h = T/N;
% x = linspace(0,T,N);
 
 

%Standard RK4
order = 4;
c = [0.0,0.5,0.5,0.0];
b = [1.0/6.0,2.0/6.0,2.0/6.0,1.0/6.0];
k = zeros(order,eqNo);
A = [0,0,0,0;
    0.5,0,0,0;
    0,0.5,0,0;
    0,0,1.0,0];

for i=2:N
        for ii=1:eqNo
            k(1,ii) = f{ii}(x(i-1),y(i-1,:));
        end
        for j = 2:order
            temp = repmat(A(j,1:j-1),eqNo,1);
            for ii=1:eqNo
                k(j,ii) = f{ii}(x(i-1) + c(j)*h, y(i-1,:) + h*dot(temp',k(1:j-1,:),1));
            end
        end
        for ii=1:eqNo
            y(i,ii) = y(i-1,ii) + h*dot(b,k(:,ii));
        end
%     k(1,1) = f{1}(x(i-1),y(i-1,:));
%     k(1,2) = f{2}(x(i-1),y(i-1,:));
%     k(2,1) = f{1}(x(i-1) + c(2)*h, [y(i-1,1)+h*A(2,1)*k(1,1), y(i-1,2)+h*A(2,1)*k(1,2)]);
%     k(2,2) = f{2}(x(i-1) + c(2)*h, [y(i-1,1)+h*A(2,1)*k(1,1), y(i-1,2)+h*A(2,1)*k(1,2)]);
%     k(3,1) = f{1}(x(i-1) + c(3)*h, [y(i-1,1)+h*(A(3,2)*k(2,1)), y(i-1,2)+h*(A(3,2)*k(2,2))]);
%     k(3,2) = f{2}(x(i-1) + c(3)*h, [y(i-1,1)+h*(A(3,2)*k(2,1)), y(i-1,2)+h*(A(3,2)*k(2,2))]);
%     k(4,1) = f{1}(x(i-1) + c(4)*h, [y(i-1,1)+h*(A(4,3)*k(3,1)), y(i-1,2)+h*(A(4,3)*k(3,2))]);
%     k(4,2) = f{2}(x(i-1) + c(4)*h, [y(i-1,1)+h*(A(4,3)*k(3,1)), y(i-1,2)+h*(A(4,3)*k(3,2))]);
%     
%     y(i,1) = y(i-1,1) + h*dot(b,k(:,1));
%     y(i,2) = y(i-1,2) + h*dot(b,k(:,2));
end
%Exact solution
t = linspace(0,4,100);
exact = y(1,1)*exp(-t).*cos(9.95*t);
plot(x,y(:,1),'r',t,exact,'o');