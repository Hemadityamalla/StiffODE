function dydt = vdp1000(t,y)
%VDP1000  Evaluate the van der Pol ODEs for mu = 1000.

dydt = [y(2); 1000*(1-y(1)^2)*y(2)-y(1)];