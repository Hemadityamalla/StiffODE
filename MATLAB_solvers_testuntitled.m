clear;
clc;

T = 3000;

opts = odeset('RelTol',1e-12,'Abstol',1e-14,'Stats','on');
fprintf("Stats for ode23s.............\n");
[t,exact] = ode23s(@vdp1000,[0 T],[2; 0],opts);%Check for the value of mu inside @vdp
plot(t,exact(:,1),'-');
hold on;
% opts = odeset('RelTol',1e-2,'Stats','on');%Using a 1 order BDF
% fprintf("\n Stats for ode15s............ \n");
% [t,exact] = ode15s(@vdp1000,[0 T],[2;0],opts);
% plot(t,exact(:,1),'-*');
% title('Solution of van der Pol Equation, \mu = 1000');
% xlabel('Time t');
% ylabel('Solution y_1');
% legend('ode23s', 'ode15s');

csvwrite('VDP_Exact.dat',[t,exact]);

