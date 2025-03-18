function u = vector_f_time(t,theta,N) 
dx = 1.5/(N+1);

e1= zeros(N,1);
e1(end)= (1+sin(5*t))/(dx) ;

u = repmat(e1,N-1,1);

e2 = zeros(N,1);
e2(1:end-1) = (1+sin(10*t))/(dx);
e2(end) =  (1+sin(5*t))/(dx) + (1+sin(10*t))/(dx);

u = [u ; e2];
