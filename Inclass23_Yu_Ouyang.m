%In this directory, you will find a derivaitve function representing a
%three gene circuit in which each gene product represses the transcription of another
%gene in a cycle. This simple three gene circuit displays oscillations and was one
%of the first synthetic circuits ever constructed (see Elowitz & Leibler
%Nature 2000). 

% 1. Start with initial conditions x(1) = 1, x(2) = x(3) = 0. 
% Run the model for 200 time units and plot the result. verify that it
% does indeed oscillate.

figure;
sol3 = ode23(@rhs3,[0 200],[1,0,0]);
plot(sol3.x,sol3.y(1,:),'r.-','LineWidth',1,'MarkerSize',18); hold on;
plot(sol3.x,sol3.y(2,:),'g.-','LineWidth',1,'MarkerSize',18); hold on;
plot(sol3.x,sol3.y(3,:),'b.-','LineWidth',1,'MarkerSize',18);
legend({'x_1','x_2','x_3'},'FontSize',18);
hold off;




% 2. What happens if you start with initial conditions x(1) = x(2) = x(3) =
% 0.5? Why?

figure;
sol3 = ode23(@rhs3,[0 200],[0.5, 0.5, 0.5]);
plot(sol3.x,sol3.y(1,:),'r.-','LineWidth',1,'MarkerSize',18); hold on;
plot(sol3.x,sol3.y(2,:),'g.-','LineWidth',1,'MarkerSize',18); hold on;
plot(sol3.x,sol3.y(3,:),'b.-','LineWidth',1,'MarkerSize',18);
legend({'x_1','x_2','x_3'},'FontSize',18);
hold off;


% There would be no oscillation, because it is at the stable point.

% 3. The strength of the repression is governed by the constant k2 which is
% the same for all the genes. What happens when you make k2 larger or
% smaller? Find the approximate value of k2 for which the system no longer
% oscillates. 
% new oscillation function to take k2 into consideration as an input


%%
%figure;
%k2 = 0.2;
%rhss = @(t,x) k2.*rhs3;
%sol3sm = ode23(@rhss,[0 200],[1,0,0]);
%plot(sol3sm.x,sol3sm.y(1,:),'r.-','LineWidth',1,'MarkerSize',18); hold on;
%plot(sol3sm.x,sol3sm.y(2,:),'g.-','LineWidth',1,'MarkerSize',18); hold on;
%plot(sol3sm.x,sol3sm.y(3,:),'b.-','LineWidth',1,'MarkerSize',18); hold on;
%k2 = 20;
%rhsl = @(t,x) k2.*rhs3;
%sol3lg = ode23(@rhsl,[0 200],[1,0,0]);
%plot(sol3lg.x,sol3lg.y(1,:),'k.-','LineWidth',1,'MarkerSize',18); hold on;
%plot(sol3lg.x,sol3lg.y(2,:),'y.-','LineWidth',1,'MarkerSize',18); hold on;
%plot(sol3lg.x,sol3lg.y(3,:),'m.-','LineWidth',1,'MarkerSize',18);
%legend({'x_1(small k2)','x_2(small k2)','x_3(small k2)',...
%    'x_1(large k2)','x_2(large k2)','x_3(large k2)'},'FontSize',18);
%hold off;

% Comment: the larger k2 is, the quiker they would reach to the stable
% point. However, the function haddle seems does not work on mine. It shows
% '"rhss" was previously used as a variable, conflicting with its use here
% as the name of a function of command'. However, when I check the class of
% rhss, it shows function_handle.

figure;
sol3sm = ode23(@rhs3ks,[0 10],[1,0,0]);
plot(sol3sm.x,sol3sm.y(1,:),'r.-','LineWidth',1,'MarkerSize',18); hold on;
plot(sol3sm.x,sol3sm.y(2,:),'g.-','LineWidth',1,'MarkerSize',18); hold on;
plot(sol3sm.x,sol3sm.y(3,:),'b.-','LineWidth',1,'MarkerSize',18); hold on;
sol3lg = ode23(@rhs3kl,[0 10],[1,0,0]);
plot(sol3lg.x,sol3lg.y(1,:),'k.-','LineWidth',1,'MarkerSize',18); hold on;
plot(sol3lg.x,sol3lg.y(2,:),'y.-','LineWidth',1,'MarkerSize',18); hold on;
plot(sol3lg.x,sol3lg.y(3,:),'m.-','LineWidth',1,'MarkerSize',18);
legend({'x_1(small k2)','x_2(small k2)','x_3(small k2)',...
    'x_1(large k2)','x_2(large k2)','x_3(large k2)'},'FontSize',18);
hold off;

% When k2 gets larger, the oscillation gets faster. 


% 4. What happens when you make k2 larger or smaller for only one of the
% genes? 

figure;
sol3sm = ode23(@rhs2s,[0 10],[1,0,0]);
plot(sol3sm.x,sol3sm.y(1,:),'r.-','LineWidth',1,'MarkerSize',18); hold on;
plot(sol3sm.x,sol3sm.y(2,:),'g.-','LineWidth',1,'MarkerSize',18); hold on;
plot(sol3sm.x,sol3sm.y(3,:),'b.-','LineWidth',1,'MarkerSize',18); hold on;
sol3lg = ode23(@rhs2l,[0 10],[1,0,0]);
plot(sol3lg.x,sol3lg.y(1,:),'k.-','LineWidth',1,'MarkerSize',18); hold on;
plot(sol3lg.x,sol3lg.y(2,:),'y.-','LineWidth',1,'MarkerSize',18); hold on;
plot(sol3lg.x,sol3lg.y(3,:),'m.-','LineWidth',1,'MarkerSize',18);
legend({'x_1(small k2 for x2)','x_2(small k2 for x2)','x_3(small k2 for x2)',...
    'x_1(large k2 for x2)','x_2(large k2 for x2)','x_3(large k2 for x2)'},'FontSize',18);
hold off;

% Yu Ouyang Comment: When k2 is small, the amplitude of the oscillation is
% larger, with a longer period. I set k2 as the changing coefficient of x2.
% The maximum of x2 is larger than those of the rest. The absolute value of
% x2 maximum equals to those of x1 and x2 minimum. In other words, the
% minimum of x2 would also be larger than those of x1 and x3, since they
% are negative. 


%%
function dx = rhs3(t,x)
dx = zeros(3,1);
dx(1) = x(2)-x(3);
dx(2) = x(3)-x(1);
dx(3) = x(1)-x(2);
end

function dx = rhs3ks(t,x,k2)
k2 = 0.2;
dx = k2.*rhs3(t,x);
end

function dx = rhs3kl(t,x,k2)
k2 = 20;
dx = k2.*rhs3(t,x);
end

function dx = rhs2s(t,x)
k2 = 0.2
dx = zeros(3,1);
dx(1) = k2.*x(2)-x(3);
dx(2) = x(3)-x(1);
dx(3) = x(1)-k2.*x(2);
end

function dx = rhs2l(t,x)
k2 = 20
dx = zeros(3,1);
dx(1) = k2.*x(2)-x(3);
dx(2) = x(3)-x(1);
dx(3) = x(1)-k2.*x(2);
end
