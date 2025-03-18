function A= matrixA(ax,bx,nx,dim)

switch dim
    case 1
        h=(bx-ax)/(nx+1);
        e=ones(nx,1);
        A=spdiags([e/h^2 (-(2/h^2))*e e/h^2],-1:1,nx,nx);
        I = speye(size(A));
        A = kron(I,A) + kron(A,I);
        J = speye(size(A));
        A = kron(J,A) + kron(A,J);
        A_dig = 0.25*diag(A);
        A_dig_mat = diag(A_dig);
        A = A - A_dig_mat;
        A = A(1:nx^3,1:nx^3);
%     case 2
        
  
end


% lap1d one dimensional finite difference approximation
% A=lap1d(n) computes a (full) matrix finite difference
% approximation of the one dimensional operator (Delta) on the
% domain Omega=(0,1) using n interior points






