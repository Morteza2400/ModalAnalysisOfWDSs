function [P_P,N_h,C_PRV,C_pump,C_TCV,C_Nor] = conversion(Pipe,Node)
%% N_h matrix for conversion
N_h=zeros(Node.N);
k=0;
KK=Node.N_tank+Node.N_joint;
for i=1:Node.N
    % if (ismember(i,Node.tank_rowAtinput))
    %     k=k+1;
    %     N_h(k,i)=1;
    % end
    if (ismember(i,Node.joint_rowAtinput))
        k=k+1;
        N_h(k,i)=1;
    end
    if (ismember(i,Node.res_rowAtinput))
        KK=KK+1;
        N_h(KK,i)=1;
    end
end

%% P_P matrix for conversion
P_P=zeros(sum(Pipe.reach));
K=0;
k=sum(Pipe.reach)-Pipe.ctr;
KK=sum(Pipe.reach)-Pipe.ctr+Pipe.Pump_n;
kk=sum(Pipe.reach)-Pipe.ctr+Pipe.Pump_n+Pipe.PRV_n;
C_PRV(1:sum(Pipe.reach),1)=0;
C_Nor(1:sum(Pipe.reach),1)=0;
C_pump(1:sum(Pipe.reach),1)=0;
C_TCV(1:sum(Pipe.reach),1)=0;
for i=1:Pipe.N
    if (ismember(i,Pipe.Normal_rowAtinput))
        for j=1:Pipe.reach(i)
        K=K+1;
        P_P(K,sum(Pipe.reach(1:i-1))+j)=1;
        C_Nor(sum(Pipe.reach(1:i-1))+j)=1;
        end
        
    end
    if (ismember(i,Pipe.pump_row))
        k=k+1;
        P_P(k,sum(Pipe.reach(1:i-1))+1)=1;
        C_pump(i)=1;
    end
    if (ismember(i,Pipe.PRV_row))
        KK=KK+1;
        P_P(KK,sum(Pipe.reach(1:i-1))+1)=1;
        C_PRV(i)=1;
    end
    if (ismember(i,Pipe.TCV_row))
        kk=kk+1;
        P_P(kk,sum(Pipe.reach(1:i-1))+1)=1;
        C_TCV(i)=1;
    end
end
end