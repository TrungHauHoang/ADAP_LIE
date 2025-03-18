function A= matrixA_nonuni(ax,bx,nx,dim)

switch dim
    case 1
        h=(bx-ax)/(nx+1);
        e=ones(nx,1);
        A=spdiags([e/h^2 (-(2/h^2))*e e/h^2],-1:1,nx-1,nx-1);
        size(A)
        temp = floor(nx/2);
        A(temp,temp) = -1/h^2;
        A(temp,temp+1) = (1/3)/h^2;
        A(temp-1,temp) = (2/3)/h^2; 
        A(temp+1,temp+1) = -1/h^2;
        A(temp+1,temp+2) = (2/3)/h^2;
        A(temp,temp+1) = (1/3)/h^2;
        
        I = speye(size(A));
        A = kron(I,A) + kron(A,I);
        size(A)
        J = speye(size(A));
        A = kron(J,A) + kron(A,J);
        A_dig = 0.25*diag(A);
        A_dig_mat = diag(A_dig);
        A = A - A_dig_mat;
        A = A(1:(nx-1)^3,1:(nx-1)^3);
%     case 2
        
  
end


% lap1d one dimensional finite difference approximation
% A=lap1d(n) computes a (full) matrix finite difference
% approximation of the one dimensional operator (Delta) on the
% domain Omega=(0,1) using n interior points






