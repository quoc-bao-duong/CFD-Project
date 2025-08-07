function [phi_L phi_T]=BCs(delta_x,n)
%{
    BCs: Boundary Conditions

    % INPUT
    delta_x,n

    % OUTPUT
    phi_L : phi LEFT
    phi_T : phi TOP
%}

    phi_L=zeros(n,1);
    phi_T=zeros(1,n);

    y=(1-delta_x/2): -delta_x : 0;

    for i=1:n
        phi_L(i,1)=1-y(i);
    end

end