function initial = function_U0(x,y,cases) % Defines initial conditions.

switch cases
    case 1
%         initial = sin((2/3)*pi*(x+0.5))*sin((2/3)*pi*(y+0.5));
        initial = sin(pi*x)*sin(pi*y);
        
    case 2
%         initial = sin((2/3)*pi*(x+0.5))*sin((2/3)*pi*(y+0.5));
    case 3
        initial = 1 + (sin((2/3)*pi*(x+0.5)))^2*(sin((2/3)*pi*(y+0.5)))^2;
%           initial = 2 + (sin(pi*x)).^2.*(sin(pi*y)).^2;
    case 4
        initial = 1 + (sin(pi*x))^2*(sin(pi*y))^2;
    case 5 
        initial = zeros(size(x));
        
end
