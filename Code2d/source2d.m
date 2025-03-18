function u = source2d(t,x,y,a)
u = ((1 + 2*a*(4/9)*pi^2)*sin((2/3)*pi*(x+0.5)).* sin((2/3)*pi*(y+0.5)))*exp(t)...
- ((2/3)*pi*(1./abs(x) + 1./abs(y)).*(sin((2/3)*pi*(x+0.5)).*cos((2/3)*pi*(y+0.5)) + cos((2/3)*pi*(x+0.5)).*sin((2/3)*pi*(y+0.5))))*exp(t);
% size(u)
% u = u(:); 

% u = exp(t)*sin((2*pi*(x + 1/2))/3).*sin((2*pi*(y + 1/2))/3) + (8*a*pi^2*exp(t)*sin((2*pi*(x + 1/2))/3).*sin((2*pi*(y + 1/2))/3))/9 - (exp(t)*(2*pi*cos((pi*(4*x + 4*y + 1))/6).*(abs(x) + abs(y))))./(3*abs(x).*abs(y));
