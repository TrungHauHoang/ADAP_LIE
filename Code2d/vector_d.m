function u = vector_d(x,y,N) 

dx = 1.5/(N+1);
% dx
e1 = zeros(N,1);
temp = y(N:N:N*N-N);
e1 = ( 1/abs(x(end))+1./abs(temp))/(dx) ;
% e1
% e1 
u = zeros(N*N-N,1);
u(N:N:N*N-N) = e1;

e2 = zeros(N,1);
e2(1:end-1) = ( 1./abs(x(1:N-1)) + 1/abs(y(end)))/(dx);
e2(end) =  (2/abs(y(end)))/(dx) + (2/abs(x(end)))/(dx);

u = [u ; e2];
