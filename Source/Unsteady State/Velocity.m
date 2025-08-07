function [u v]=Velocity(delta_x,n)
%{
INPUT
   delta_x,n
OUTPUT
    u: Velocity (x-axis)
    v: Velocity (y-axis)
%}

    u=zeros(n); % matrix u with size: n x n
    v=zeros(n); % matrix v with size: n x n

    x= delta_x/2 :delta_x : 1;
    y= (1-delta_x/2) : -delta_x : 0
    for i=1:n
        u(:,i)=x(i);
        v(i,:)=-y(i);
    end

end