function f = dfdt(t,y,k)

kJac = dkdT(temperature(t));

L = 2.5*y(1) + 16*y(3) + y(2) + y(4) + y(5) + y(6) + y(7) + y(8);
M = 2.5*y(1) + 6*y(3) + y(2) + y(4) + y(5) + y(6) + y(7) + y(8);
dk4dT = (kJac(4,1)*k(4,2))/(k(4,2) + k(4,1)*L) - (L*k(4,1)*(kJac(4,1)*k(4,2) - kJac(4,2)*k(4,1)))/(k(4,2) + k(4,1)*L)^2;
dk12dT = (kJac(12,2)*k(12,1)*M)/(k(12,2) + k(12,1)*M) + ...
         (M*k(12,2)*(kJac(12,1)*k(12,2) - kJac(12,2)*k(12,1)))/(k(12,2) + k(12,1)*M)^2;

f(1) = kJac(2,2) - kJac(2,1) + kJac(3,2) - kJac(3,1) + kJac(9,1) - kJac(9,2) - kJac(11,1);
f(2) = kJac(1,2) - kJac(1,1) - dk4dT + kJac(6,1) - kJac(6,2) + kJac(7,1) + kJac(10,1);
f(3) = kJac(3,1) - kJac(3,2) + kJac(7,1) + kJac(8,1) - kJac(8,2); 
f(4) = kJac(1,2) - kJac(1,1) + kJac(2,1) - kJac(2,2) + kJac(3,1) - kJac(3,2) - dk4dT - kJac(5,1) - kJac(6,1) + kJac(6,2) + kJac(11,1);
f(5) = kJac(1,1) - kJac(1,2) - kJac(2,2) + kJac(2,2);
f(6) = kJac(1,2) - kJac(1,2) + kJac(2,1) - kJac(2,2) - kJac(3,1) + kJac(3,2) + 2*kJac(5,1) - kJac(7,1) - kJac(8,1) + kJac(8,2) + 2*dk12dT;
f(7) = dk4dT - kJac(5,1) - kJac(6,1) + kJac(6,2) - kJac(7,1) - 2*kJac(10,1) - kJac(11,1);
f(8) = kJac(10,1) + kJac(11,1);


f = f.*dTdt(t);
%f = zeros(size(y));

end