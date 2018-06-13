clear;
clc;
h = 0.2;
tfinal = 4.0;
n = int64(tfinal/h);

y = zeros(n,1);
t = linspace(0,tfinal,n);

y(1,1) = 1;

%Starting step for Explicit Adams-Bashfort(2)
y(2,1) = y(1,1) + h*dydt(t(1),y(1,1));

for i=1:n-2
    y(i+2) = y(i+1) + 0.5*h*(3*dydt(t(i+1),y(i+1,1)) - dydt(t(i),y(i)));
    
end
plot(t,y);





function f = dydt(t,y)

f = (1 - 2*t)*y;

end