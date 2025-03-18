function u = vector_e(x,y,z,N,theta)

dx = 1.5/(N+1);
e1 = sparse(N,1);
u = sparse(N*N*N,1);
temp = y(N:N:N*N-N);
temp1 = x(1:N);

for i = 1:N
    u_temp = sparse(N*N-N,1);
    e1 = ( 1/abs(x(end)) + 1./abs(temp(1:N-1)) + 1/abs(temp1(i)) - theta(i) )/(dx) ;
    u_temp(N:N:N*N-N) = e1;
    e2 = sparse(N,1);
    e2(1:end-1) = ( 1./abs(x(1:N-1)) + 1/abs(y(end)) + 1/abs(temp1(i)) - theta(i))/(dx);
    e2(end) =  (2/abs(y(end)) + 2/abs(x(end)) + 2/abs(temp1(i)) - 2*theta(i))/(dx);
    u_temp = [u_temp ; e2];
    u(N*N*(i-1)+1:i*N*N) = u_temp;
end

e1 = ( 1./abs(x(1:N*N)) + 1./abs(y(1:N*N)) + 1/abs(temp1(end)) - theta(end) )/(dx) ;
u(N*N*N-N*N+1:N*N*N) = u(N*N*N-N*N+1:N*N*N) + e1;


