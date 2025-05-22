function plotGraphWithCustomLabels(Node, incidenceMatrix, nodeValues, edgeValues, ModeNumber, real_p, imag_p)
    % Create the directed graph from the incidence matrix
    [numNodes, numEdges] = size(incidenceMatrix);
    edges = zeros(numEdges, 2);
    figure
    % Build the edge list from the incidence matrix
    for e = 1:numEdges
        nodes = find(incidenceMatrix(:, e));
        edges(e, 1) = nodes(incidenceMatrix(nodes, e) == 1);
        edges(e, 2) = nodes(incidenceMatrix(nodes, e) == -1);
    end
    
    G = digraph(edges(:,1), edges(:,2));
    
    % Get edge order for labeling
    edgeTable = G.Edges;
    edgeMatrix = table2array(edgeTable);
    [~, order] = ismember(edgeMatrix, edges, 'rows');
    
    % Custom node numbering as 'N1', 'N2', ..., without numeric labels
    nodeLabels = arrayfun(@(x) ['N' num2str(x)], 1:numNodes, 'UniformOutput', false);
    
    % Custom edge numbering for labeling
    edgeLabels = arrayfun(@(x) ['E' num2str(x)], order, 'UniformOutput', false); 

    % Plot the graph
    nodePositions = [Node.xdata', Node.ydata']; % Node positions
    
    % Plot the graph without numeric node labels
    p = plot(G, 'XData', nodePositions(:,1), 'YData', nodePositions(:,2), ...
        'MarkerSize', 3, 'LineWidth', 2, 'NodeLabel', []); % Disable numeric node labels
    
    % Assign node and edge values to the plot
    p.NodeCData = nodeValues;
    p.EdgeCData = edgeValues(order); % Order the edge values correctly
    p.NodeColor = 'flat';
    p.EdgeColor = 'flat';
    
    % Customize edge labels appearance
    % p.EdgeLabel = edgeLabels;    % Set edge labels
    % p.EdgeLabelColor = 'black';  % Set edge label color to black
    
    % Set colormap and colorbar
    colormap('jet');
    colorbar;
    
    % Customize the title with mode number, real and imaginary parts
    title(['PF for mode with ', ' , real = ', num2str(real_p), ' , imag = ', num2str(imag_p), ' (', num2str(imag_p/2/pi),' Hz)']);
    
    % Label the axes
    xlabel('Node position X (m)');
    ylabel('Node position Y (m)');
    
    % Adjust axis properties
    axis tight; 
    axis equal; 
    axis manual;

    % Add custom node labels manually (placing the labels near the nodes)
    % labelOffset = 0.1;  % Adjust this value to control the label offset
    % 
    % for i = 1:numNodes
    %     % Only display the custom node labels (N1, N2, ...) and no numeric labels
    %     text(nodePositions(i,1) - labelOffset, nodePositions(i,2) - labelOffset, ...
    %         nodeLabels{i}, 'FontSize', 7, 'Color', 'black', ...
    %         'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
    % end
end
