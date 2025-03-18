function u = vector_f(theta,N)

dx = 1.5/(N+1);
% dx
e1 = sparse(N,1);
u = sparse(N*N*N,1);
for i = 1:N
    u_temp = sparse(N*N-N,1);
    e1 = theta(i)/dx;
    u_temp(N:N:N*N-N) = e1;
    e2 = sparse(N,1);
    e2(1:end-1) = theta(i)/(dx);
    e2(end) =  (2*(theta(i)))/(dx);
    u_temp = [u_temp ; e2];
    u(N*N*(i-1)+1:i*N*N) = u_temp;
end

e1 = theta(end)/(dx) ;
u(N*N*N-N*N+1:N*N*N) = u(N*N*N-N*N+1:N*N*N) + e1;


