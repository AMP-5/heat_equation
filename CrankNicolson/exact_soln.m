function [u] = exact_soln(x, t)
    %DATA Summary of this function goes here
    %   Detailed explanation goes here
    Nt = length(t);
    Nx = length(x);
    u= zeros(Nx,Nt);
    for p=1:Nt
        for q=1:Nx-1
            u(q,p) =exp(-t(p))*sin(x(q));
        end    
    end 
    
end

