function u = vector_e(x,y,N,theta) 
dx = 1.5/(N+1);
% dx
e1 = zeros(N,1);
temp = y(N:N:N*N-N);
e1 = ( 1/abs(x(end)) + 1./abs(temp) - 1)/(dx) ;
% e1
% e1 
u = zeros(N*N-N,1);
u(N:N:N*N-N) = e1;


e2 = zeros(N,1);
e2(1:end-1) = (1*( 1./abs(x(1:N-1))+1/abs(y(end)) - 1))/(dx);
e2(end) =  (2/abs(y(end)))/(dx) + (2/abs(x(end)))/(dx) - 2/dx;

u = [u ; e2];

