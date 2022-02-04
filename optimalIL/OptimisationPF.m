%Optimisation method Particle formation
fig=figure;
set(fig,'Position',[501,165,904,715]);
set(gcf,'color','w');
ax2=axes();%subplot(2,1,2);
hold(ax2,'on');
xlabel(ax2,'$t$','Interpreter','Latex','Fontsize',14);
ylabel(ax2,'$x_i(t)$','Interpreter','Latex','Fontsize',14);

rng('default');
rng(5);

Nrobots=2;
n=2*Nrobots;
m=0.2*ones(1,Nrobots);
gamma=ones(1,Nrobots);
d=ones(1,Nrobots-1);
D=zeros(n,n); %noise width of the langevin equations

x0=zeros(1,Nrobots*2);
%Initial conditions
for i=1:2:2*Nrobots
  D(i+1,i+1)=0.5;  
  x0(i)=(1/Nrobots)*i+((1/Nrobots)*(i+1)-(1/Nrobots)*i)*rand; 
end
tspan=[0 10];

R=ga(@(gains) ParticleControlSimulation(gains,m,gamma,d,D,x0,Nrobots,n,tspan),(Nrobots-1)*2,[],[],[],[],0,5);

rpos=zeros(length(t),Nrobots);
for k=1:length(t)
    x=mvnrnd(M(k,:),reshape(yx(k,:),n,n));
    x(2:2:end)=[];
    rpos(k,:)=x;
end
M(:,2:2:end)=[];
for i=1:Nrobots
    %plot(ax1,t,M(:,i),'k:')
    plot(ax2,t,rpos(:,i),'k:')
end