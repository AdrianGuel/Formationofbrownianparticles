%Model based on Guel-Cortez,2021
%https://link.springer.com/chapter/10.1007/978-3-030-82064-0_16
function [A,B,U]=GenerateSyS(m,kp,kd,gamma,d)
    A=zeros(length(m)*2,length(m)*2);
    B=zeros(length(m)*2,length(m)*2);
    U=zeros(length(m)*2,1);
    %set A and B ones
    g=1;    
    for i=1:length(m)*2
        for k=i+1:length(m)*2
            if (mod(k,2)==0 && mod(i,2)~=0)
                A(i,k)=1;
                B(i+1,k)=1;
                U(k)=1;
                break;
            end
        end
    end
    for i=2:2:length(m)*2
        s1=1;
        s2=1;
        if g<=length(kp)
            for k=i-1:length(m)*2
                if s1==3
                    break;
                elseif mod(k,2)~=0
                    A(i,k)=((-1)^s1)*kp(g)/m(g);
                    A(i,k+1)=((-1)^s1)*kd(g)/m(g);
                    s1=s1+1;
                    if B(i,k+1)==1
                        B(i,k+1)=-kp(g)/m(g);
                        U(k+1)=d(g);
                    end
                end
            end
            for k=i-1:length(m)*2-1
                if s2==2
                    break;
                elseif mod(k,2)~=0
                    A(i,k+1)=A(i,k+1)-gamma(g)/m(g);
                    s2=s2+1;
                end
            end            
        end
            g=g+1;
    end
    A(end,end)=-gamma(end)/m(end);
end