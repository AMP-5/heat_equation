function [U] = approx_soln_U(x, t)
    %DATA Summary of this function goes here
    %   Detailed explanation goes here
    Nt = length(t);
    Nx = length(x);
    k = (t(end)-t(1))/(Nt-1);
    h = (x(end)-x(1))/(Nx-1);
    U= zeros(Nx,Nt);
    b = ones(Nx-2,1);
    
    for q=1:Nx-1
            U(q,1) = sin(x(q));
    end
    L1 = spdiags([(-b/(2*h^2)) (b/k+b/(h^2)) (-b/(2*h^2))], -1:1, Nx-2, Nx-2); % construction of sparse matrix L
    
    L2 = spdiags([(b/(2*h^2)) (b/k-b/(h^2))  (b/(2*h^2))], -1:1, Nx-2, Nx-2); % construction of sparse matrix L
    
    for n=2:Nt
    
    U(2:Nx-1,n) = L1\(L2*U(2:Nx-1,n-1));
 
    end
end


