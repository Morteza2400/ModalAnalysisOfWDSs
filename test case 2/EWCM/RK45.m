    function [Y] = RK45(t,X,f,h)   %RK45
        K1=f(t,X);
        K2=f(t+h/2,X+h/2*K1);
        K3=f(t+h/2,X+h/2*K2);
        K4=f(t+h,X+h*K3);
        Y=X+h/6*(K1+2*K2+2*K3+K4);
    end