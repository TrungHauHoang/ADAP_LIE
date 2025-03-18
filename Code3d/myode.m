function dydt = myode(t,y,A,C,choice,X,Y,Z,a,nx,theta,normvector,Cy,Cz)
switch choice
    case 1
        dydt = A*y + C*y + a*vector_c(nx) + vector_d(X,Y,Z,nx);
    case 2
        dydt = C*y + vector_d(X,Y,Z,nx);
    case 3
        dydt = A*y + a*vector_c(nx);
    case 4
        dydt = C*y + vector_f(theta,nx);
    case 5
        dydt = A*y + a*vector_c(nx) + C*y + vector_e(X,Y,Z,nx,theta);
    case 6
        dydt = A*y + C*y + Cy*y + Cz*y + (2-1./normvector).*y - 1./normvector;
    case 7
        dydt = C*y + Cy*y + Cz*y ;
    case 8
        dydt = A*y + (2-1./normvector).*y - 1./normvector ;  
    case 9
        dydt = C*y + Cy*y + Cz*y ;
    case 10
        dydt = A*y + C*y + Cy*y + Cz*y + (2-1./normvector).*y - 1./normvector ;
        
end
end
