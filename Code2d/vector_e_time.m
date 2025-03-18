function u = vector_e_time(t,x,y,N,theta) 
dx = 1.5/(N+1);
% dx
e1 = zeros(N,1);
temp = y(N:N:N*N-N);
e1 = ((1+sin(5*t))*( 1/abs(x(end)) + 1./abs(temp) - 1))/(dx) ;
% e1
% e1 
u = zeros(N*N-N,1);
u(N:N:N*N-N) = e1;


e2 = zeros(N,1);
e2(1:end-1) = ((1+sin(10*t))*( 1./abs(x(1:N-1))+1/abs(y(end)) - 1))/(dx);
e2(end) =  ((1/abs(x(end)) + 1/abs(y(end))-1)*(+2+sin(5*t)+sin(10*t)))/(dx);

u = [u ; e2];

