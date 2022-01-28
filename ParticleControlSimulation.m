%Stochastic simulation of line formation of brownian particles
%Adrian Guel Cortez 2022
%https://link.springer.com/chapter/10.1007/978-3-030-82064-0_16

fig=figure;
set(fig,'Position',[501,165,904,715]);
set(gcf,'color','w');
% ax1=subplot(2,1,1);
% hold(ax1,'on');
% xlabel(ax1,'$t$','Interpreter','Latex','Fontsize',14);
% ylabel(ax1,'$\mu(t)$','Interpreter','Latex','Fontsize',14);

ax2=axes();%subplot(2,1,2);
hold(ax2,'on');
xlabel(ax2,'$t$','Interpreter','Latex','Fontsize',14);
ylabel(ax2,'$x_i(t)$','Interpreter','Latex','Fontsize',14);

rng('default');
rng(5);

Nrobots=12;
m=0.2*rand(1,Nrobots);
gamma=4*rand(1,Nrobots);
d=ones(1,Nrobots-1);
kp=5*rand(1,Nrobots-1);
kd=5*rand(1,Nrobots-1);
[A,B,u]=GenerateSyS(m,kp,kd,gamma,d);
u(end)=0; %no force in the leading particle.
n=2*Nrobots;
D=zeros(n,n); %noise width of the langevin equations

x0=zeros(1,Nrobots*2);
%Initial conditions
for i=1:2:2*Nrobots
  D(i+1,i+1)=0.5;  
  x0(i)=(1/Nrobots)*i+((1/Nrobots)*(i+1)-(1/Nrobots)*i)*rand; 
end


tspan=[0 10];
opts=odeset('RelTol',1e-8,'AbsTol',1e-10);
%Mean Value
[t,M]=ode45(@(t,y) Mux(t,y,A,B,u),tspan,x0,opts);
%Covariance Matrix
[t,yx]=ode45(@(t,y) Sigmax(t,y,A,D,n),t,D,opts);

%Stochastic simulation from normal distirbution
rpos=zeros(length(t),Nrobots);
for k=1:length(t)
    x=mvnrnd(M(k,:),reshape(yx(k,:),n,n));
    x(2:2:end)=[];
    rpos(k,:)=x;
end
%M(:,2:2:end)=[];
for i=1:Nrobots
    %plot(ax1,t,M(:,i),'k:')
    plot(ax2,t,rpos(:,i),'k:')
end

function dydt=Mux(t,y,A,B,u)
    dydt=A*y+B*u;
end

function dydt=Sigmax(t,y,A,D,n)
    At=A';
    aux=reshape(y,n,n);
    dydt=kron(eye(n,n),A)*aux(:)+kron(eye(n,n),aux)*At(:)+2*D(:);
end