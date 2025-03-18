function u = vector_f(theta,N) 
dx = 1.5/(N+1);

e1= zeros(N,1);
e1(end)= 1/(dx) ;

u = repmat(e1,N-1,1);

e2 = zeros(N,1);
e2(1:end-1) = 1/(dx);
e2(end) =  1/(dx) + 1/(dx);

u = [u ; e2];
