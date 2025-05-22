function [f_q]=friction(Pipe,g,q,fv,Fluid,R)
if (Fluid.Friction==0)
Re = abs(q)./(Fluid.kv*pi*(Pipe.DD')/4)+.01;
B1 = (-2.457.*log((7./Re).^0.9+0.27*(Pipe.ee')./(Pipe.DD'))).^16;
B2 = (37530./Re).^16;
fl = 8*((8./Re).^12+(B1+B2).^(-3/2)).^(1/12);
else
    fl =Pipe.ee';
end
Fw = fl.*R';
Fv = fv/2/g./(Pipe.AA').^2;
Fm = zeros (1,sum(Pipe.reach))';
Sj = ones(1,sum(Pipe.reach))';
f_q = Sj .* (Fw + Fm + Fv);   % replace q
end