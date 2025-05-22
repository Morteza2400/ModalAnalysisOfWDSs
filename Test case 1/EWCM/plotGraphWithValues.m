function plotGraphWithValues(Node, incidenceMatrix, nodeValues, edgeValues, ModeNumber, real_p, imag_p,i)
    % Create the directed graph from the incidence matrix
    [numNodes, numEdges] = size(incidenceMatrix);
    edges = zeros(numEdges, 2);
    
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

     % Customize edge labels appearance

    % p.EdgeLabelColor = 'black';  % Set edge label color to black
    
    % Plot the graph
    nodePositions = [Node.xdata', Node.ydata']; % Node positions
    
    % Plot graph with formatted edges and nodes
    p = plot(G, 'XData', nodePositions(:,1), 'YData', nodePositions(:,2), ...
        'MarkerSize', 6, 'LineWidth', 2, 'NodeLabel', []); % Disable numeric node labels
    
    % Assign node and edge values to the plot
    p.NodeCData = nodeValues;
    p.EdgeCData = edgeValues(order); % Order the edge values correctly
    p.NodeColor = 'flat';
    p.EdgeColor = 'flat';
    
    % Apply colormap (consistent across subplots)
    colormap('jet');
    p.EdgeLabel = edgeLabels;    % Set edge labels
    p.EdgeLabelColor = 'black';  % Set edge label color to black
    p.EdgeFontSize = 9;         % Set edge label font size
    % Customize title with mode information
    title(['PF for mode', ' ', num2str(i)], ...
           'FontSize', 12, 'FontWeight', 'Bold', 'FontName', 'Times New Roman');

    % Set axes labels
    xlabel('Node position X (m)', 'FontSize', 12, 'FontName', 'Times New Roman');
    % ylabel('Node position Y (m)', 'FontSize', 12, 'FontName', 'Times New Roman');
    
    % Adjust axis properties
    axis tight; 
    axis equal; 
    axis manual;

    % Add custom node labels manually with offsets to prevent overlap
    labelOffset = 0.1;  % Adjust this value to control label placement
    
    for i = 1:numNodes
        text(nodePositions(i,1) - labelOffset, nodePositions(i,2) - labelOffset, ...
            nodeLabels{i}, 'FontSize', 9, 'FontName', 'Times New Roman', 'FontWeight', 'Bold', ...
            'Color', 'black', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
    end
end
