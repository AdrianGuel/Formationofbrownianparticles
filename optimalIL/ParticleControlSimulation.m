%Stochastic simulation of line formation of brownian particles
%Adrian Guel Cortez 2022
%https://link.springer.com/chapter/10.1007/978-3-030-82064-0_16
function IL=ParticleControlSimulation(gains,m,gamma,d,D,x0,Nrobots,n,tspan)

    kp=gains(1:Nrobots-1);
    kd=gains(Nrobots:end);

    [A,B,u]=GenerateSyS(m,kp,kd,gamma,d);
    u(end)=0; %no force in the leading particle.

    opts=odeset('RelTol',1e-8,'AbsTol',1e-10);
    %Mean Value
    [t,M]=ode45(@(t,y) Mux(t,y,A,B,u),tspan,x0,opts);
    %Covariance Matrix
    [t,yx]=ode45(@(t,y) Sigmax(t,y,A,D,n),t,D,opts);    
    
    E=zeros(1,length(t));
    for k=1:length(t)
        E(k)=InformationEnergy(A,M(k,:)',B,u,reshape(yx(k,:),n,n),D);
    end
    
    IL=trapz(t,sqrt(E));
    
    function Et=InformationEnergy(A,Mu,B,u,Sigma,D)
        Et=(A*Mu+B*u)'*(Sigma\(A*Mu+B*u))+0.5*trace((Sigma\(A*Sigma+Sigma*A'+2*D))^2);
    end

    function dydt=Mux(t,y,A,B,u)
        dydt=A*y+B*u;
    end

    function dydt=Sigmax(t,y,A,D,n)
        At=A';
        aux=reshape(y,n,n);
        dydt=kron(eye(n,n),A)*aux(:)+kron(eye(n,n),aux)*At(:)+2*D(:);
    end
end