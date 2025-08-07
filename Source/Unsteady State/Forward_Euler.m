function [Phi_new]=Forward_Euler(n,rho,delta_x,delta_y,Gamma,u,v,phi_T,phi_L,delta_t,t_max,Phi_initial)

K1=Gamma*delta_y/delta_x;
K2=Gamma*delta_x/delta_y;
K3=rho*delta_y;
K4=rho*delta_x;


Phi_old=Phi_initial;
Phi_new=zeros(n,n);

for t=delta_t:delta_t:t_max
           % Cell (1,1)
           A_1=-(3*K1+3*K2)*Phi_old(1,1)-(0.5*K3*u(1,1)-K1)*Phi_old(1,2)-(-0.5*K4*v(1,1)-K2)*Phi_old(2,1);
           A_2=2*K1*phi_L(1,1)+2*K2*phi_T(1,1)+0.5*K3*u(1,1)*phi_L(1,1)-0.5*K4*v(1,1)*phi_T(1,1);
           Phi_new(1,1)=Phi_old(1,1)+(A_1+A_2)*delta_t/(rho*delta_x*delta_y);

           %Cell (1,n)
           B_1=-(-0.5*K3*u(1,n)-K1)*Phi_old(1,n-1)-(0.5*K3*u(1,n)+K1+3*K2-0.5*K4*v(1,n))*Phi_old(1,n)-(-0.5*K4*v(1,n)-K2)*Phi_old(2,n);
           Phi_new(1,n)=Phi_old(1,n)+(B_1)*delta_t/(rho*delta_x*delta_y);

           % Cell (n,1)
           C_1=-(0.5*K4*v(n,1)-K2)*Phi_old(n-1,1)-(3*K1+K2-0.5*K4*v(n,1))*Phi_old(n,1)-(-K1+0.5*K3*v(n,1))*Phi_old(n,2);
           C_2=0.5*K3*u(n,1)*phi_L(n,1)+2*K1*phi_L(n,1);
           Phi_new(n,1)=Phi_old(n,1)+(C_1+C_2)*delta_t/(rho*delta_x*delta_y);

           % Cell (n,n)
           D_1=-(-0.5*K3*u(n,n)-K1)*Phi_old(n,n-1)-(0.5*K3*u(n,n)+K1+K2-0.5*K4*v(n,n))*Phi_old(n,n)-(K4*v(n,n)-K2)*Phi_old(n-1,n);
           Phi_new(n,n)=Phi_old(n,n)+(D_1)*delta_t/(rho*delta_x*delta_y);

           % Cell top (1,j)
           for j=2:n-1
           E_1=-(-0.5*K3*u(1,j)-K1)*Phi_old(1,j-1)-(2*K1+3*K2)*Phi_old(1,j)-(-0.5*K4*v(1,j)-K2)*Phi_old(2,j)-(0.5*K3*u(1,j)-K1)*Phi_old(1,j+1)-0.5*K4*v(1,j)*phi_T(1,j)+2*K2*phi_T(1,j);
           Phi_new(1,j)=Phi_old(1,j)+(E_1)*delta_t/(rho*delta_x*delta_y);
           end

           % Cell bottom (n,j)
           for j=2:n-1
           F_1=-(-0.5*K3*u(n,j)-K1)*Phi_old(n,j-1)-(2*K1+K2-0.5*K4*v(n,j))*Phi_old(n,j)-(0.5*K4*v(n,j)-K2)*Phi_old(n-1,j)-(0.5*K3*u(n,j)-K1)*Phi_old(n,j+1);
           Phi_new(n,j)=Phi_old(n,j)+(F_1)*delta_t/(rho*delta_x*delta_y);
           end

           % Cell left (i,1)
           for i=2:n-1
           G_1=-(0.5*K4*v(i,1)-K2)*Phi_old(i-1,1)-(3*K1+2*K2)*Phi_old(i,1)-(0.5*K3*u(i,1)-K1)*Phi_old(i,2)-(-0.5*K4*v(i,1)-K2)*Phi_old(i+1,1)+2*K1*phi_L(i,1)+0.5*K3*u(i,1)*phi_L(i,1);
           Phi_new(i,1)=Phi_old(i,1)+(G_1)*delta_t/(rho*delta_x*delta_y);
           end

           % Cell right (i,n)
           for i=2:n-1
           M_1=-(-0.5*K3*u(i,n)-K1)*Phi_old(i,n-1)-(0.5*K4*v(i,n)-K2)*Phi_old(i-1,n)-(0.5*K3*u(i,n)+K1+2*K2)*Phi_old(i,n)-(-0.5*K4*v(i,n)-K2)*Phi_old(i+1,n);
           Phi_new(i,n)=Phi_old(i,n)+(M_1)*delta_t/(rho*delta_x*delta_y);
           end 

           % Cell center (i,j)
           for i=2:n-1
               for j=2:n-1
               N_1=-(-0.5*K3*u(i,j)-K1)*Phi_old(i,j-1)-(-K2+0.5*K4*v(i,j))*Phi_old(i-1,j)-(2*K1+2*K2)*Phi_old(i,j)-(0.5*K3*u(i,j)-K1)*Phi_old(i,j+1)-(-0.5*K4*v(i,j)-K2)*Phi_old(i+1,j);
               Phi_new(i,j)=Phi_old(i,j)+(N_1)*delta_t/(rho*delta_x*delta_y);
               end
           end

           Phi_old=Phi_new;
end
end
