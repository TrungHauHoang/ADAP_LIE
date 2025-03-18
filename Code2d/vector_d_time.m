function u = vector_d_time(t,x,y,N) 

b1 = 1+sin(5*t);
b2 = 1+sin(10*t);

dx = 1.5/(N+1);
e1 = zeros(N,1);
temp = y(N:N:N*N-N);
e1 = (b1*(1/abs(x(end))+1./abs(temp)))/(dx) ;

u = zeros(N*N-N,1);
u(N:N:N*N-N) = e1;

e2 = zeros(N,1);
e2(1:end-1) = (b2*(1./abs(x(1:N-1))+1/abs(y(end))))/(dx);
e2(end) =  ((1/abs(x(end)) + 1/abs(y(end)))*(2+sin(5*t)+sin(10*t)))/(dx);

u = [u ; e2];






% function u = vector_d_time(t,x,y,N) 
% 
% dx = 1.5/(N+1);
% e1 = zeros(N,1);
% temp = y(N:N:N*N-N);
% e1 = ((1+sin(5*t))*(1/abs(x(end))+1./abs(temp)))/(dx) ;
% 
% u = zeros(N*N-N,1);
% u(N:N:N*N-N) = e1;
% 
% e2 = zeros(N,1);
% e2(1:end-1) = ((1+sin(10*t))*(1./abs(x(1:N-1))+1/abs(y(end))))/(dx);
% e2(end) =  ((1/abs(x(end)) + 1/abs(y(end)))*(+2+sin(5*t)+sin(10*t)))/(dx);
% 
% u = [u ; e2];
