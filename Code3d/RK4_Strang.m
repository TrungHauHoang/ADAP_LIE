function [ w1 ] = RK4_Strang(t,x,w0)
    
    w1 = zeros(length(w0),1);
    
           
    for k = 1:length(t)-1
        dt = t(k+1)-t(k);
        dw1 = w0 - w0.^3;
        dw2 = (w0+(dt/2)*dw1) - (w0+(dt/2)*dw1).^3;
        dw3 = (w0+(dt/2)*dw2) - (w0+(dt/2)*dw2).^3;
        dw4 = (w0+dt*dw3) - (w0+dt*dw3).^3;
        w1 = w0+(dt/6)*(dw1+2*dw2+2*dw3+dw4);
        w0 = w1;
    end

end