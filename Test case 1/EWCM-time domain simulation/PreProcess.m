function [Pipe,Node,Fluid,A,A1,A2,Res_row,B,B1,q0,h0,h1,...
    h2,q,L,R,C,M_p,M_h,P_P,N_h,C_PRV,C_pump,C_TCV,C_Nor,RR,A_save,N_h_total]=PreProcess(Pipe,Fluid,Node,g)
C=0;
M_p = sum(Pipe.reach);
M_h=Node.N+sum(Pipe.reach-1);
A = zeros (M_h,M_p);
x1 = 1; x2 = 1;
for i = 1: Pipe.N
    A (Pipe.up(i),x1) = 1;
    A (Pipe.down(i),x1 + Pipe.reach(i) - 1 ) = -1;
    A (Node.N+x2:Node.N+x2+Pipe.reach(i)-2,x1:x1+Pipe.reach(i)-1) = fX (Pipe.reach(i));
    x1 = x1 + Pipe.reach(i);
    x2 = x2 + Pipe.reach(i) -1;
end

x2 = 1;
for j=1:Pipe.N
    if (Pipe.reach(j)>1)
        for k=1:Pipe.reach(j)-1
            Node.xdata(Node.N+k+x2-1)=Node.xdata(Pipe.up(j))+k*(Node.xdata(Pipe.down(j))-Node.xdata(Pipe.up(j)))/Pipe.reach(j);
            Node.ydata(Node.N+k+x2-1)=Node.ydata(Pipe.up(j))+k*(Node.ydata(Pipe.down(j))-Node.ydata(Pipe.up(j)))/Pipe.reach(j);
        end
    end
    x2 = x2 + Pipe.reach(i) -1;
end

A_save=A;
x1 = 1;
for i = 1:Pipe.N
    for j = x1:x1+Pipe.reach(i)-1   % miyad L C R ro baraye har Pipe hesab mikone va mesh ro dar nazar migire
        L(j,1) = Pipe.l(i)./Pipe.reach(i)./g./Pipe.A(i);
        C(j) = g*Pipe.A(i)*Pipe.l(i)/Pipe.reach(i)/Pipe.a(i)^2;
        R(j) = Pipe.l(i)/Pipe.reach(i)/2/g./Pipe.D(i)/Pipe.A(i)^2;
        RR(j) = Pipe.e(i);
        Pipe.LL(j,1) = Pipe.l(i);
    end
    x1 = x1 + Pipe.reach(i);
end

Res_row = zeros (Node.N_Res,sum(Pipe.reach));      % Shift end node row to the end of A
% to in for miyad aval satr haee ke tosh reseivoir dare ro save mikone bad
% onsatr haro az A hazfmikone va mifrestateshon akhare A
for i = 1 : Node.N_Res
    Res_row (i,:) = A (Node.res_rowAtinput(i)-(i-1),:);
    A (Node.res_rowAtinput(i)-(i-1),:) = [];
    A (end+1,:) = Res_row (i,:);
end

B = abs(A);

A1 = A(1:M_h-Node.N_Res,:);
A2 = A(M_h-Node.N_Res+1:end,:);
B1 = B(1:M_h-Node.N_Res,:);
Pipe.ctr=Pipe.Pump_n+Pipe.PRV_n+Pipe.TCV_n;
q0 = zeros (M_p,1);
h0 = zeros (M_h,1);
[P_P,N_h,C_PRV,C_pump,C_TCV,C_Nor,N_h_total] = conversion(Pipe,Node,A,A_save);
h0(1:Node.N) = Node.H0';

Node.Q0_demand=N_h*Node.Q0_dem';
Node.d=N_h*Node.d';
x1 = 1; x2 = Node.N;
for i =1 : Pipe.N
    q0(x1:x1+Pipe.reach(i)-1) = Pipe.ini_q (i);
    Pipe.DD(x1:x1+Pipe.reach(i)-1) = Pipe.D(i);
    Pipe.ee(x1:x1+Pipe.reach(i)-1) = Pipe.e(i);
    Pipe.AA(x1:x1+Pipe.reach(i)-1) = Pipe.A(i);
    for j = 1:Pipe.reach(i)-1
        h0 (x2+j) = Node.H0(Pipe.up(i)) - j* (Node.H0(Pipe.up(i)) - Node.H0(Pipe.down(i)))/Pipe.reach(i);
    end
    x1 = x1 + Pipe.reach(i);
    x2 = x2 + Pipe.reach(i) -1 ;
end
q = q0;
h1 = h0;
h2 = zeros (Node.N_Res,1);
for i = 1: Node.N_Res % inja ham miyad satr reservoir ro mibare akhar
    h1(Node.res_rowAtinput(i)-(i-1)) = [];
    h2(i) = Node.H0(Node.res_rowAtinput(i));
end
end