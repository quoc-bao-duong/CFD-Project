function [A,B]=Central(n,rho,delta_x,delta_y,Gamma,u,v,phi_T,phi_L)

% Create Matrix A size n^2 x n^2
A=zeros(n*n);
B=zeros(n*n,1);
% 
K1=Gamma*delta_y/delta_x;
K2=Gamma*delta_x/delta_y;
K3=rho*delta_y;
K4=rho*delta_x;

% Cell Phi(1,1) (Equation 1)
A(1,1)=3*K1+3*K2+K3*u(1,1);
A(1,2)= 0.5*K3*u(1,1)-K1;
A(1,n+1)=-0.5*K4*v(1,1)-K2;
B(1,1)=2*K1*phi_L(1,1)+2*K2*phi_T(1,1)+0.5*K3*u(1,1)*phi_L(1,1)-0.5*K4*v(1,1)*phi_T(1,1);

% Cell Phi(1,N) (Equation N)
A(n,n-1)=-0.5*K3*u(1,n)-K1;
A(n,n)=0.5*K3*u(1,n)+K1+3*K2-0.5*K4*v(1,n);
A(n,2*n)=-0.5*K4*v(1,n)-K2;
B(n,1)=0;

% Cell Phi(N,1) (Equation (N-1)*N+1)
A((n-1)*n+1,(n-2)*n+1)=0.5*K4*v(n,1)-K2
A((n-1)*n+1,(n-1)*n+1)=K3*u(n,1)+3*K1+K2-0.5*K4*v(n,1);
A((n-1)*n+1,(n-1)*n+2)=-K1+0.5*K3*v(n,1);
B((n-1)*n+1,1)=0.5*K3*u(n,1)*phi_L(n,1)+2*K1*phi_L(n,1);

% Cell Phi(N,N) (Equation N*N)
A(n*n,n*n-1)=-0.5*K3*u(n,n)-K1;
A(n*n,n*n)=0.5*K3*u(n,n)+K1+K2-0.5*K4*v(n,n);
A(n*n,n*n-n)=K4*v(n,n)-K2;
B(n*n,1)=0;

% Cell TOP Phi
for j=2:n-1
    A(j,j)=2*K1+3*K2;
    A(j,j-1)=-0.5*K3*u(1,j)-K1;
    A(j,j+1)=0.5*K3*u(1,j)-K1;
    A(j,j+n)=-0.5*K4*v(1,j)-K2;
    B(j,1)=-0.5*K4*v(1,j)*phi_T(1,j)+2*K2*phi_T(1,j);
end

% Cell BOT Phi
for j=2:n-1
    A((n-1)*n+j,(n-1)*n+j)=2*K1+K2-0.5*K4*v(n,j);
    A((n-1)*n+j,(n-1)*n+j-1)=-0.5*K3*u(n,j)-K1;
    A((n-1)*n+j,(n-1)*n+j+1)=0.5*K3*u(n,j)-K1;
    A((n-1)*n+j,(n-1)*n+j-n)=-K2+0.5*K4*v(n,j);
    B((n-1)*n+j,1)=0;
end

% Cell LEFT Phi
for i=2:n-1
    A((i-1)*n+1,(i-1)*n+1)=3*K1+2*K2;
    A((i-1)*n+1,(i-1)*n+2)=0.5*K3*u(i,1)-K1;
    A((i-1)*n+1,(i-2)*n+1)=0.5*K4*v(i,1)-K2;
    A((i-1)*n+1,i*n+1)=-0.5*K4*v(i,1)-K2;
    B((i-1)*n+1,1)=2*K1*phi_L(i,1)+0.5*K3*u(i,1)*phi_L(i,1);
end

% Cell RIGHT Phi
for i=2:n-1
    A(i*n,i*n)=0.5*K3*u(i,n)+K1+2*K2;
    A(i*n,i*n-1)=-0.5*K3*u(i,n)-K1;
    A(i*n,i*n-n)=0.5*K4*v(i,n)-K2;
    A(i*n,i*n+n)=-0.5*K4*v(i,n)-K2;
    B(i*n,1)=0;
end

% Inner Cell Phi
for i=2:n-1
    for j=2:n-1
        A((j-1)*n+i,(j-1)*n+i)=2*K1+2*K2;
        A((j-1)*n+i,(j-1)*n+i-1)=-0.5*K3*u(i,j)-K1;
        A((j-1)*n+i,(j-1)*n+i+1)=0.5*K3*u(i,j)-K1;
        A((j-1)*n+i,(j-1)*n+i-n)=-K2+0.5*K4*v(i,j);
        A((j-1)*n+i,(j-1)*n+i+n)=-0.5*K4*v(i,j)-K2;
    end
end


end