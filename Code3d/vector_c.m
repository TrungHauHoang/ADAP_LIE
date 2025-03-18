function u = vector_c(N)
dx = 1.5/(N+1);

c2 = sparse(N,1);
c2(1) = 2/(dx^2);
c2(end) = 2/(dx^2);
c2(2:end-1) = 1/(dx^2);

c3 = sparse(N,1);
c3(1) = 1/(dx^2);
c3(end) = 1/(dx^2);

u = repmat(c3,N-2,1);
u = [c2 ; u ; c2];
u = repmat(u,N-2,1);

d2 = sparse(N,1);
d2(1) = 3/(dx^2);
d2(end) = 3/(dx^2);
d2(2:end-1) = 2/(dx^2);

temp = repmat(c2,N-2,1);
temp = [d2;temp;d2];
u = [temp ; u ; temp];

