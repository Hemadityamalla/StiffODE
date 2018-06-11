clear;
clc;

TOL = 1e-4;
ATOL = 1e-2;

N = 25;
t = linspace(0,4,N+1);
h = (4-0)/N;

y = zeros(N+1,2);

y(1,1) = 1;
y(1,2) = 1;

for i=1:N-1
    %Computing the first term using RK2
    if i == 1
        k1 = f(t(i),y(i,:));
        ytemp = y(i,:) + 0.5*h*k1;
        k2 = f(t(i)+0.5*h,ytemp);
        
        y(i+1,:) = y(i,:) + h*k2; 
    end
    
    %Solving the system using Newton-Raphson
    
    %Initial guess using RK-2
    
    k1 = f(t(i+1),y(i+1,:));
    ytemp = y(i+1,:) + 0.5*h*k1;
    k2 = f(t(i+1)+0.5*h,ytemp);
    yGuess = y(i+1,:) + h*k2;
    
    F = (4/3)*y(i+1,:) - (1/3)*y(i,:) + (2*h/3)*f(t(i+2),yGuess) - yGuess;
    J = [(2*h/3)-1, -(2*h/3)-1; (4*h/3)-1, -2*h-1];
    
    yGuess2 = [0,0];
    iter = 1;
    while max(abs(F)) > ATOL && (max(abs(yGuess - yGuess2)) - TOL) > ATOL
        if iter ~= 1
            yGuess = yGuess2;
        end
        yGuess2 = yGuess - (J\F')';
        iter = iter+1;
    end
    
    y(i+2,:) = yGuess2;
    
end

plot(t,y(:,1),t,(t+1).*exp(-t));















function func = f(t,y)

func = [y(1) - y(2), 4*y(1) - 3*y(2)];

end