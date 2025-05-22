    function [Xm] = fX(x)   % Xm
        Xm = zeros (x-1,x);
        for j = 1: x-1
            Xm(j,j) = -1;
            Xm(j,j+1)=1;
        end
    end