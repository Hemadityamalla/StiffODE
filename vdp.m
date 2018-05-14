function dydt = vdp(t,y)
%VDP  Evaluate the van der Pol ODEs for a given mu
mu = 3.0;
dydt = [y(2); mu*(1-y(1)^2)*y(2)-y(1)];