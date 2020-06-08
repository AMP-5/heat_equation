function [U] = approx_soln_U(x, t)
    %DATA Summary of this function goes here
    %   Detailed explanation goes here
    Nt = length(t);
    Nx = length(x);
    k = (t(end)-t(1))/(Nt-1);
    h = (x(end)-x(1))/(Nx-1);
    U= zeros(Nx,Nt);
    
    for q=1:Nx-1
            U(q,1) = sin(x(q));
    end
    U(2:Nx-1,2) = (1-(2*k/(h^2)))*U(2:Nx-1,1)+(k/(h^2))*(U(3:Nx,1)+U(1:Nx-2,1));            
    for n=3:Nt
       U(2:Nx-1,n) = U(2:Nx-1,n-2)-(4*k/(h^2))*U(2:Nx-1,n-1)+(2*k/(h^2))*(U(3:Nx,n-1)+U(1:Nx-2,n-1));            
    end
end