function dydt = vdp1000(t,y)
%VDP1000  Evaluate the scaled van der Pol ODEs for mu = 1e-6.
 
dydt = [y(2); 1e6*((1-y(1)^2)*y(2)-y(1))];