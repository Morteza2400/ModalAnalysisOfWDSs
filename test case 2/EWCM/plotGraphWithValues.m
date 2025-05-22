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
    
    % Get edge order for labeling
    edgeTable = G.Edges;
    edgeMatrix = table2array(edgeTable);
    [~, order] = ismember(edgeMatrix, edges, 'rows');
    
    % Custom labels
    nodeLabels = arrayfun(@(x) ['N' num2str(x)], 1:numNodes, 'UniformOutput', false);
    edgeLabels = arrayfun(@(x) ['E' num2str(x)], order, 'UniformOutput', false);

    % Create a figure with white background and specific size
    fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 7, 4.5]);  % [left bottom width height]

    % Node positions
    nodePositions = [Node.xdata', Node.ydata'];

    % Plot the graph
    p = plot(G, 'XData', nodePositions(:,1), 'YData', nodePositions(:,2), ...
        'MarkerSize', 6, 'LineWidth', 2, 'NodeLabel', []);
    p.NodeCData = nodeValues;
    p.EdgeCData = edgeValues(order);
    p.NodeColor = 'flat';
    p.EdgeColor = 'flat';

    % Customize appearance
    colormap('jet');
    c = colorbar;
    c.Label.String = 'PF magnitude';
    c.Label.FontSize = 11;

    set(gca, 'XTick', [], 'YTick', []);
    set(gca, 'XColor', 'none', 'YColor', 'none');
    axis tight; 
    axis equal; 
    axis manual;

    % Add custom node labels
    labelOffset = 0.7;
    for i = 1:numNodes
        text(nodePositions(i,1) - labelOffset, nodePositions(i,2) - labelOffset, ...
            nodeLabels{i}, 'FontSize', 8, 'Color', 'black', ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
    end
set(gca, 'FontName', 'Times New Roman');
set(findall(gcf,'-property','FontName'), 'FontName', 'Times New Roman');
    % Save as high-resolution BMP with mode number
filenameBMP = sprintf('Graph_Mode%d.bmp', ModeNumber);
print(fig, filenameBMP, '-dbmp', '-r600');  % Save as BMP at 600 DPI



filenamePNG = sprintf('Graph_Mode%d.png', ModeNumber);
exportgraphics(fig, filenamePNG, 'Resolution', 600, 'BackgroundColor', 'white');
    close(fig);  % Optional: close figure if you're doing this in a loop
end
