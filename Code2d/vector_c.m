function u = vector_c(N)
dx = 1.5/(N+1);

c2= zeros(N,1);
c2(1) = 1/(dx^2);
c2(end)= 1/(dx^2);
u =repmat(c2,N-2,1);

d2 = zeros(N,1);
d2(1) = 2/(dx^2);
d2(end) = 2/(dx^2);
d2(2:end-1) = 1/(dx^2);

u = [d2 ; u ; d2];
