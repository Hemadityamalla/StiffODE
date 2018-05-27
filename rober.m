function dydt = rober(t,y)
%robertson chemical kinetics problem.

dydt = [(-0.04*y(1) + 1e4*y(2)*(3)); (0.04*y(1) - 1e4*y(2)*y(3) - 3e7*y(2)^2); (3e7*y(2)^2)];
