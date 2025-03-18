function u = vector_c_time(t,N)

b1 = 1+sin(5*t);
b2 = 1+sin(10*t);

dx = 1.5/(N+1);

c2= zeros(N,1);
c2(1) = (b1)/(dx^2);
c2(end)= (b1)/(dx^2);
u =repmat(c2,N-2,1);

d2 = zeros(N,1);
d2(1) = (b1 + b2)/(dx^2);
d2(end) = (b1 + b2)/(dx^2);
d2(2:end-1) = b2/(dx^2);

u = [d2 ; u ; d2];
