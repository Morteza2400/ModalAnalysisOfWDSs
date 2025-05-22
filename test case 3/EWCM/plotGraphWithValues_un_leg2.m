function plotGraphWithCustomLabels2(Node, incidenceMatrix, nodeValues, edgeValues, ModeNumber, real_p, imag_p)
    % Create the directed graph from the incidence matrix
    [numNodes, numEdges] = size(incidenceMatrix);
    edges = zeros(numEdges, 2);
    % Define figure size (adjust as needed)
figWidth = 15;  % Width in cm
figHeight = 16; % Height in cm

% Create figure with specific size
figure('Units', 'centimeters', 'Position', [5, 5, figWidth, figHeight]);

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
    
    % Set colormap and colorbar
    colormap('jet');
    colorbar;
    xlabel(''); % Remove X-axis label
    ylabel(''); % Remove Y-axis label
    set(gca, 'XTick', [], 'YTick', []); % Remove tick marks
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Remove axis numbers
    axis off; % Completely remove the axis
    % Customize the title with mode number, real and imaginary parts
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
    title(['Weighted Average PF for Modes Below fc']);
    set(gcf, 'Color', 'w')
    set(gca, 'Color', 'w')
    % Label the axes
    % xlabel('Node position X (m)');
    % ylabel('Node position Y (m)');
    
    % Adjust axis properties
    axis tight; 
    axis equal; 
    axis manual;

    %% Increase Size of Specific Nodes
    specificNodes = [1,4,489,367, 414, 450, 179]; % Nodes to make larger
    highlight(p, specificNodes, 'MarkerSize', 6); % Increase size (default was 3)

specificNodes2 = [1];  % Nodes you want to mark with a square

% Get X and Y positions of the specific nodes
xSpecial = nodePositions(specificNodes2, 1);
ySpecial = nodePositions(specificNodes2, 2);

% Highlight size (optional)
highlight(p, specificNodes2, 'MarkerSize', 6);

% Overlay a square marker on the node
hold on
plot(xSpecial, ySpecial, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k'); % 'k' = black, 's' = square
hold off

    %% Add Custom Labels Only for Specific Nodes
    customLabels = {'Reservoir','N4','N_A','N17', 'N18', 'N19', 'N10'}; % Custom labels for each node
    labelOffset = 10.4;  % Adjust this value to control the label offset

    for i = 1:length(specificNodes)
        nodeIndex = specificNodes(i);
        text(nodePositions(nodeIndex,1) - 5*labelOffset, nodePositions(nodeIndex,2) + labelOffset, ...
            customLabels{i}, 'FontSize', 10, 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'Color', 'black', ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    end
% Define filenames
filenamePDF = 'Graph_HighQuality.pdf';
filenameBMP = 'Graph_HighQuality.jpg';

% Save as high-quality PDF (vector format, preserves figure size)
exportgraphics(gcf, filenamePDF, 'ContentType', 'vector', 'Resolution', 600, 'BackgroundColor', 'white');

% Save as high-quality BMP (raster format)
exportgraphics(gcf, filenameBMP, 'Resolution', 600, 'BackgroundColor', 'white');
end
