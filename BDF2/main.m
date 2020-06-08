clear all; close all; format short

Nr=6; %No of refinements

for r=1:Nr
    xL=0;xR=pi;T=1;
    Nx = 16*2^r;
    %Step size h
    h = (xR-xL)/Nx;
    %As h=O(k)
    Nt = ceil(T/(0.1*(h)));
    %Dividing [0,pi] into Nx eql parts to get x
    x=linspace(xL,xR,Nx+1);
    x=x(:);
    %Dividing [0,1] into Nt eql parts to get t
    t=linspace(0,T,Nt+1);
    t=t(:);
    %Step size k
    k = (t(end)-t(1))/Nt;
    [u] = exact_soln(x,t);
    [U] = approx_soln_U(x,t);
    abs_error = abs(u-U);
    Error(r) = max(max(abs_error));
    
end
Error
R = log2(Error(1:end-1)./Error(2:end))
