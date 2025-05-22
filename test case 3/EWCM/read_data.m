function [t_total,deltaT,Fluid,Node,Pipe,g,PRV,TCV,Pump]=read_data(A)
sizem=size(A);
rows_A=sizem(1);
for i=1:rows_A
    if (strcmp(A{i,1},'[Hydrulic OPTIONS]'))
        HY_OPTIONS_row = i;
    end
    if (strcmp(A{i,1},'[NODES]'))
        NODES_row=i;
    end
    if (strcmp(A{i,1},'[PIPES]'))
        PIPES_row = i;
    end
    if (strcmp(A{i,1},'[END]'))
        END_row = i;
    end
end
deltaT=cell2mat(A(HY_OPTIONS_row+1,2)); %number of time step
g=cell2mat (A(HY_OPTIONS_row+2,2));
t_total=cell2mat(A(HY_OPTIONS_row+3,2)); %Temperature of pipes network (K)
Fluid.kv=cell2mat(A(HY_OPTIONS_row+4,2));
Fluid.c=cell2mat(A(HY_OPTIONS_row+5,2));
Fluid.density=cell2mat(A(HY_OPTIONS_row+6,2));
Fluid.Friction=cell2mat(A(HY_OPTIONS_row+7,2));
%------------------------
% Number of Nodes
j = 0;
for i= NODES_row + 4:PIPES_row - 1
    if (A{i,2}>-10^20)
        j = j+1;
    end
end
Node.N = j; %Number of Nodes
fprintf('The number of Nodes is %d\n', Node.N); %Get Number of Pipes
j = 0;
for i=PIPES_row+4:END_row-1
    if (A{i,4} > 0)
        j = j+1;
    end
end
Pipe.N =j; %Number of pipes
fprintf('The number of pipes is %d\n', Pipe.N);
j=0;
%% reading data from Node section
KK=0;
NN=0;
MM=0;
Node.res_rowAtinput=0;
Node.tank_rowAtinput=0;
Node.joint_rowAtinput=0;
Node.Tank_save=[];
for i=NODES_row+4:NODES_row+4+Node.N-1
    j = j+1;
    Node.ID(j)=j;
    S1 = num2str (A{i, 1});
    Node.string_ID{j} = num2str(A{i,1});
    if (strcmp(A{i,3}, 'Reservoir'))
        KK=KK+1;
        Node.Res(KK)=j;
        Node.res_rowAtinput(KK)=j;
        Node.Res_save=Node.Res;
    end
    if (strcmp(A{i,3}, 'Tank'))
        MM=MM+1;
        Node.Tank(MM)=j;
        Node.Tank_save=Node.Tank;
        Node.tank_rowAtinput(MM)=j;
        Node.At(MM)=cell2mat(A(i,6));
        Node.Tank_ini_h(MM)=cell2mat(A(i,7));
        Node.Tank_demand(MM)=cell2mat(A(i,4));
    end
    if (strcmp(A{i,3}, 'Joint'))
        NN=NN+1;
        Node.Nc(NN)=j;
        Node.joint(NN)=j;
        Node.joint_rowAtinput(NN)=j;
        Node.joint_save=Node.joint;
    end
    Node.type{j} = S1;
    Node.Elevation(j)=cell2mat(A(i,2));
    Node.Q0_dem(j)=cell2mat(A(i,4));
    Node.H0(j)=cell2mat(A(i,5));
    Node.xdata(j)=cell2mat(A(i,6));
    Node.ydata(j)=cell2mat(A(i,7));
    Node.newID{j} = num2str(A{i,1});

end
Node.N_Res=KK;
Node.N_tank=MM;
Node.N_joint=NN;
Node.Ntot=(1:Node.N);
Node.H0_save=Node.H0;
if (Node.N_tank>0)
    Node.Tank_save=Node.Tank;
end

Node.Res_save=Node.Res;
%% reading data from Pipe section
KK=0;
NN=0;
MM=0;
DD=0;
hh=0;
j = 0;
Pipe.valve_rowAtinput=0;
Pipe.pump_rowAtinput=[];
Pipe.PRV_rowAtinput=[];
Pipe.TCV_rowAtinput=[];
Pipe.pump_row=[];
Pipe.PRV_row=[];
Pipe.TCV_row=[];
Pipe.Normal_rowAtinput=0;
for i= PIPES_row + 4:PIPES_row + 4 + Pipe.N-1
    j = j+1;
    S1 = num2str(A{i,4});
    Pipe.ID(j)=j;
    Pipe.up(j)=cell2mat (A(i,2));
    Pipe.down(j)=cell2mat (A(i,3));
    Pipe.type{j} = num2str(A{i,4});
    Pipe.D(j) = cell2mat (A(i,5));
    Pipe.l(j)=cell2mat(A(i,6));
    Pipe.ini_q(j) = cell2mat(A(i,7));
    Pipe.e(j)=cell2mat(A(i,8));
    Pipe.reach(j)=cell2mat(A(i,9));
    Pipe.a(j)=cell2mat(A(i,10));
    if (strcmp(S1, 'Pump'))
        KK=KK+1;
        Pipe.pump(KK)=j;
        Pipe.pump_row(KK)=j;
        Pipe.pump_rowAtinput(KK)=sum(Pipe.reach(1:j));
        Pump.v_set(KK)=cell2mat (A(i, 11));
        Pump.T(KK)=cell2mat (A(i, 12));
        Pump.A(KK)=cell2mat (A(i, 13));
        Pump.B(KK)=cell2mat (A(i, 14));
        Pump.C(KK)=cell2mat (A(i, 15));
        Pump.S_ini(KK)=cell2mat (A(i, 16));
    end
    if (strcmp(S1, 'PRV'))
        MM=MM+1;
        Pipe.PRV(MM)=j;
        Pipe.PRV_row(MM)=j;
        Pipe.PRV_rowAtinput(MM)=sum(Pipe.reach(1:j));
        PRV.alpha_o(MM)=cell2mat (A(i, 11));
        PRV.alpha_c(MM)=cell2mat (A(i, 12));
        PRV.h_set(MM)=cell2mat (A(i, 13));
        PRV.Xm_ini(MM)=cell2mat (A(i, 14));
        PRV.Cv_A(MM)=cell2mat (A(i, 15));
        PRV.Cv_B(MM)=cell2mat (A(i, 16));
        PRV.Cv_C(MM)=cell2mat (A(i, 17));
    end
    if (strcmp(S1, 'TCV'))
        NN=NN+1;
        Pipe.TCV(NN)=j;
        Pipe.TCV_row(NN)=j;
        Pipe.TCV_rowAtinput(NN)=sum(Pipe.reach(1:j));
        TCV.Kv(NN)=cell2mat (A(i, 11));
        TCV.D(NN)=cell2mat (A(i, 12));
    end
    if (strcmp(S1, 'Normal'))
        DD=DD+1;
        Pipe.Normal_rowAtinput(DD)=j;
        Pipe.Normal_row(DD)=sum(Pipe.reach(1:j));
    end
end
Pipe.PRV_n=MM;
Pipe.TCV_n=NN;
Pipe.Pump_n=KK;
Pipe.Normal_n=DD;
Pipe.valve_n=hh;
PRV.num=MM;
Pump.num=KK;
TCV.num=NN;
%%
for y=1:Pipe.N
    Pipe.A(y) = pi*Pipe.D(y)*Pipe.D(y)/4;
    Pipe.V(y)=Pipe.A(y)*Pipe.l(y);
end

