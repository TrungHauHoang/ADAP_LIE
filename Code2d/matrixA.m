function A= matrixA(nx,dim)

switch dim
    case 1
        h=1.5/(nx+1);
        e=ones(nx,1);
        A=spdiags([e/h^2 (-(2/h^2))*e e/h^2],-1:1,nx,nx);
        I = speye(size(A));
        A = kron(I,A) + kron(A,I);
        
    case 2
        h=1/(nx+1);
        e=ones(nx,1);
        A=spdiags([e/h^2 -2/h^2*e e/h^2],-1:1,nx,nx);
        
        A(1,1) = (-2/3)/h^2;
        A(1,2) = (2/3)/h^2;
        A(end,end) = (-2/3)/h^2;
        A(end,end-1) = (2/3)/h^2;
        I = speye(size(A));
        A = kron(I,A) + kron(A,I);
end


% lap1d one dimensional finite difference approximation
% A=lap1d(n) computes a (full) matrix finite difference
% approximation of the one dimensional operator (Delta) on the
% domain Omega=(0,1) using n interior points






