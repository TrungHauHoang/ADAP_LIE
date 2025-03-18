function dydt = myode(t,y,d,A,C,choice,X,Y,a)
switch choice
    case 1
        dydt = A*y + C*y + d;
    case 2
        dydt = A*y + C*y + source2d(t,X,Y,a);
end
end

