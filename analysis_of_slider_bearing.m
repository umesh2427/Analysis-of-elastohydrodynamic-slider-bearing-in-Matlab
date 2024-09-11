N=50;
M=50;
C=0.030;
% K=0.0;
W=0.0;
T0=0;
Z1=25;
X1=10;
alpha=X1/Z1;
delxbar=1/N;
delzbar=1/M;
force =0.0;
eetastar=1;
const1=(X1^2)/(Z1^2);
ITER=1000;
for I=1:N+1
    for J=1:M+1
        p(I,J)=0.0;
        T(I,J)=T0;
    end
end


for K=1:ITER
    sumij=0.0;
    for I=2:N
        X(I)=(1/N)*(I-1);
        h(I)=(2/3)*(2-X(I));
        hm=(2/3)*(2-(X(I)-0.5*delxbar));
        hp=(2/3)*(2-(X(I)+0.5*delxbar));        
        hm1=(2/3)*(2-(X(I)-delxbar));        
        hp1=(2/3)*(2-(X(I)+delxbar));
        cubh=h(I)^3;
        cubhm=hm^3;
        cubhp=hp^3;
        const2=(cubhp+cubhm+2*const1*cubh);
        A=(const1*cubh)/const2;
        CA=cubhp/const2;
        D=cubhm/const2;
        E=(.5*delxbar)*(C/X1)*((hp1-hm1)/const2);
        for J=2:M
            Z(J)=(1/M)*(J-1);
            p(I,J)=A*p(I,J+1)+A*p(I,J-1)+CA*p(I+1,J)+D*p(I-1,J)-E;
            qx=(h(I)-cubh*((p(I,J)-p(I-1,J))/delxbar))/2;
            qy=-(cubh/2)*((p(I,J)-p(I,J-1))/delzbar);
            aa=(1/delxbar)+(qy*alpha)/(qx*delzbar);
            bb=(2*eetastar)/(qx*h(I));
            cc=(6*h(I))/(eetastar*qx);
            T(I,J)=(((qy*alpha)/(qx*delzbar))*T(I,J-1)+(1/delxbar)*T(I-1,J)+bb+cc*(((p(I,J)-p(I-1,J))/delxbar)^2+alpha^2*(((p(I,J)-p(I,J-1))/delzbar)^2)))/aa;
            sumij=sumij+p(I,J);
           
        end 
       
    end
    sum(K+1)=sumij;
    percentage = abs(sum(K+1)-sum(K))/abs(sum(K+1));
    if percentage < 0.0001
        break
    end
end

Y=K
surf(T);
cc=max(max(p))
