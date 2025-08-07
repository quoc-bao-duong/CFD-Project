function A=Result(Phi,n,Phi_T,Phi_L)
%{ 
 OUTPUT:
    
 INPUT

%}

% Convert matrix size n^2 x 1 into matrix size n x n

    Phi_new=reshape(Phi, n, n)';
    A=Phi_new;
    
end