% function plotGraphWithCustomLabels(Node, incidenceMatrix, nodeValues, edgeValues, ModeNumber, real_p, imag_p)
%     % Create the directed graph from the incidence matrix
%     [numNodes, numEdges] = size(incidenceMatrix);
%     edges = zeros(numEdges, 2);
% 
%     for e = 1:numEdges
%         nodes = find(incidenceMatrix(:, e));
%         edges(e, 1) = nodes(incidenceMatrix(nodes, e) == 1);
%         edges(e, 2) = nodes(incidenceMatrix(nodes, e) == -1);
%     end
% 
%     G = digraph(edges(:,1), edges(:,2));
% 
%     % Custom edge numbering for labeling
%     edgeLabels = arrayfun(@(x) ['E' num2str(x)], 1:numEdges, 'UniformOutput', false); 
%     nodeLabels = arrayfun(@(x) ['' num2str(x)], 1:numNodes, 'UniformOutput', false);
% 
%     % Plot the graph
%     % figure;
%     nodePositions = [Node.xdata', Node.ydata']; % Node positions
% 
%     % Plot the graph with customized labels
%     p = plot(G, 'XData', nodePositions(:,1), 'YData', nodePositions(:,2), ...
%         'NodeLabel', nodeLabels, 'MarkerSize', 6, 'LineWidth', 2);
% 
%     % Assign node and edge values to the plot
%     p.NodeCData = nodeValues;
%     p.EdgeCData = edgeValues;
%     p.NodeColor = 'flat';
%     p.EdgeColor = 'flat';
% 
%     % Customize node labels appearance (e.g., black color for node labels)
%     p.NodeLabel = nodeLabels;
%     p.NodeLabelColor = 'black';  % Set node label color to black
% 
%     % Customize edge labels appearance (e.g., red color for edge labels)
%     p.EdgeLabel = edgeLabels;    % Set the edge labels
%     p.EdgeLabelColor = 'red';    % Set edge label color to red
% 
%     % Set colormap and colorbar
%     colormap('jet');
%     colorbar;
% 
%     title(['PF - mode num = ', num2str(ModeNumber), ' , real = ', num2str(real_p), ' , imag = ', num2str(imag_p), ' (', num2str(imag_p/2/pi),' Hz)']);
%     xlabel('Node position X');
%     ylabel('Node position Y');
%     axis tight; 
%     axis equal; 
%     axis manual
% end


function plotGraphWithCustomLabels(Node, incidenceMatrix, nodeValues, edgeValues, ModeNumber, real_p, imag_p)
    % Create the directed graph from the incidence matrix
    [numNodes, numEdges] = size(incidenceMatrix);
    edges = zeros(numEdges, 2);
    
    for e = 1:numEdges
        nodes = find(incidenceMatrix(:, e));
        edges(e, 1) = nodes(incidenceMatrix(nodes, e) == 1);
        edges(e, 2) = nodes(incidenceMatrix(nodes, e) == -1);
    end
    
    G = digraph(edges(:,1), edges(:,2));
    
    % Custom edge numbering for labeling
    edgeLabels = arrayfun(@(x) ['E' num2str(x)], 1:numEdges, 'UniformOutput', false); 
    nodeLabels = arrayfun(@(x) ['' num2str(x)], 1:numNodes, 'UniformOutput', false);

    % Plot the graph
    nodePositions = [Node.xdata', Node.ydata']; % Node positions
    
    % Plot the graph with customized labels
    p = plot(G, 'XData', nodePositions(:,1), 'YData', nodePositions(:,2), ...
        'NodeLabel', nodeLabels, 'MarkerSize', 6, 'LineWidth', 2);
    
    % Assign node and edge values to the plot
    p.NodeCData = nodeValues;
    p.EdgeCData = edgeValues;
    p.NodeColor = 'flat';
    p.EdgeColor = 'flat';
    
    % Customize node labels appearance (e.g., black color for node labels)
    p.NodeLabel = nodeLabels;
    p.NodeLabelColor = 'black';  % Set node label color to black
    
    % Customize edge labels appearance (e.g., red color for edge labels)
    p.EdgeLabel = edgeLabels;    % Set the edge labels
    p.EdgeLabelColor = 'red';    % Set edge label color to red
    
    % Increase the font size of node and edge labels
    p.NodeFontSize = 9;         % Set node label font size
    p.EdgeFontSize = 9;         % Set edge label font size
    
    % Set colormap and colorbar
    colormap('jet');
    colorbar;
    
    % Customize the title with mode number, real and imaginary parts
    title(['PF - mode num = ', num2str(ModeNumber), ' , real = ', num2str(real_p), ' , imag = ', num2str(imag_p), ' (', num2str(imag_p/2/pi),' Hz)']);
    
    % Label the axes
    xlabel('Node position X');
    ylabel('Node position Y');
    
    % Adjust axis properties
    axis tight; 
    axis equal; 
    axis manual;
end
