A = [-5,6;2,-4];
b = [-3;10];

x = A\b;

fprintf("The values of u and v are %f, %f \n",x(1),x(2));

V = [2,3,4,5,7,10];
I = [5.2,7.8,10.7,13,19.3,27.5];
R = (sum(dot(V,I))/sum(dot(I,I)))









function FINALPROBLEM
T0 = 0.0;
TFINAL = 1.0;
x1_0 = 1.0;
x2_0 = 0.0;
[t,x] = ode45(@myode,[T0 TFINAL],[x1_0; x2_0]);
plot(t,x(:,1));
xlabel('Time t');
ylabel('Solution x_1');
end


function dydt = myode(t,x)

dydt = [x(2); -2*x(2) - (6*x(1)^3 + 4*x(1))];

end








