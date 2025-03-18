function dydt = myode_time(t,y,A,C,choice,X,Y,a,nx,theta)
switch choice
    case 1
        dydt = A*y + C*y + a*vector_c(nx) + vector_d(X,Y,nx);
    case 2
        dydt = C*y + vector_d(X,Y,nx);
    case 3
        dydt = A*y + a*vector_c(nx);
    case 4
        dydt = C*y + vector_f(theta,nx);
    case 5
        dydt = A*y + a*vector_c(nx) + C*y + vector_e(X,Y,nx,theta);
    case 6
        dydt = A*y + C*y + a*vector_c_time(t,nx) + vector_d_time(t,X,Y,nx);
    case 7
        dydt = C*y + vector_d_time(t,X,Y,nx);
    case 8
        dydt = A*y + a*vector_c_time(t,nx);
    case 9
        dydt = C*y + vector_f_time(t,theta,nx);
    case 10
        dydt = A*y + a*vector_c_time(t,nx) + C*y + vector_e_time(t,X,Y,nx,theta);
end
end

