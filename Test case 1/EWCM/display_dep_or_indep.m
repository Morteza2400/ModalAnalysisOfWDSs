function display_dep_or_indep(dep,Pipe)

figure(1)
directed_G = digraph(Pipe.up,Pipe.down);
edges(:,1)=Pipe.up;
edges(:,2)=Pipe.down;
% Extract the edge table
edgeTable = directed_G.Edges;

% Convert the edge table to a matrix
edgeMatrix = table2array(edgeTable);
% Find the order of appearance of rows in the first matrix
[~, order] = ismember(edgeMatrix, edges, 'rows');

for i=1:Pipe.N
    if ismember(order(i),dep)
        id(i)=0;
    else
        id(i)=1;
    end
end
plot(directed_G, 'EdgeLabel',id)
title '0 = Dependent    1 = Independent'
end 