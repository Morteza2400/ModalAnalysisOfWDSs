function [P_P,N_h,C_PRV,C_pump,C_TCV,C_Nor,N_h_total] = conversion(Pipe,Node,A,A_save)
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
N_h_total=find_permutation_matrix(A,A_save);
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
function P = find_permutation_matrix(A, B)
    % This function finds the permutation matrix P that transforms A into B
    % where A and B are adjacency matrices of a graph.
    % A: The original adjacency matrix.
    % B: The permuted adjacency matrix.
    % P: The permutation matrix such that B = P * A.
    
    % Get the size of the matrices
    n = size(A, 1);
    
    % Initialize the permutation matrix with zeros
    P = zeros(n);
    
    % Loop through each row of B to find the corresponding row in A
    for i = 1:n
        for j = 1:n
            if isequal(A(i, :), B(j, :))
                P(j, i) = 1;
                break;
            end
        end
    end
end
end