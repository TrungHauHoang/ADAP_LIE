function A=lap3d_uy(ax,bx,internalPoints,cases)

% A=lap1d(n) computes a (full) matrix mid point
% approximation of the one dimensional operator (\nabla) on the
% domain Omega=(0,1) using n interior points, 

if cases == 1 % Homogeneous Dirichlet boundary condition
    dof = internalPoints; 
    h = (bx-ax)/(dof+1);    
    A_nabla = lap1d_nabla(ax,bx,dof,cases);
    Iz = speye(dof);
    A = kron(Iz,A_nabla);
    A = kron(A,Iz);
elseif cases == 2 % Periodic boundary condition

    dof = internalPoints+1; 
    A_nabla = lap1d_nabla(dof,cases);
    Iz = speye(dof);
    A = (kron(Iz,A_nabla) + kron(A_nabla,Iz));
    
end
