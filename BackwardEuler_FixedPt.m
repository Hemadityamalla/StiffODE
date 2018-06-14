clear;
clc;

h = 0.2;
tfinal = 5.0;
n = int64(tfinal/h);

y = zeros(n,1);
t = linspace(0,tfinal,n);

%Initial Condition
y(1,1) = 0.2;

for i=1:n-1
    %Solve for y(i+1) using fixedPt iteration
    initGuess = y(i,1);
    parameters = [h, y(i,1),t(i+1)]';
    sol = fixedPt(initGuess, parameters);
    y(i+1,1) = sol;
end
plot(t,y);

function func = f(t,y)

    func = 2*y*(1-y);
    
end

function root = fixedPt(guess, params)
    h = params(1,1);
    yn = params(2,1);
    t = params(3,1);
    TOL = 1e-6;
    MAXITER = int64(1e6);
    iter = 0;
    root = 0;
    while abs(root - guess) > TOL && iter < MAXITER
       if iter ~= 0
           guess = root;
       end
       root = yn + h*f(t,guess);
       iter = iter + 1;
    end
    if iter > MAXITER
        fprintf("No solution found!\n");
    end
end
