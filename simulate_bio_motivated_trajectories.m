num_genes = 20000;
num_jumps = 20;
module_size_range = [10, 200];


% initialize the tree of life
way_points = zeros(num_genes,1);
terminal_nodes = 1;
adj = 1;
directed_edges = [];
branching_prob = 0.2;


% simulating the tree
for i=1:num_jumps
    % weights for selecting branching point
    if i==1
        weights = 1;
    else
        G = digraph(directed_edges(:,2), directed_edges(:,1));
        d = distances(G);
        dist_to_terminal = sum(d(terminal_nodes,:),1);
        weights = zeros(1,size(adj,1));
        weights(terminal_nodes) = (1-branching_prob)/length(terminal_nodes);
        weights(setdiff(1:end,terminal_nodes)) = exp(-dist_to_terminal(setdiff(1:end,terminal_nodes)));
        weights(setdiff(1:end,terminal_nodes)) = weights(setdiff(1:end,terminal_nodes))/sum(weights(setdiff(1:end,terminal_nodes)))*branching_prob;
    end
    
    % pick a point for branching
    branching_node = randsample(1:size(way_points,2), 1, true, weights);
    adj(branching_node,end+1)=1;
    adj(end+1,branching_node)=0;
    directed_edges = [directed_edges; [branching_node, size(adj,1)]];
    
    % pick a subset of genes within the size range
    randomed_size = randsample(min(module_size_range):max(module_size_range),1);
    genes_involved = randsample(1:num_genes, randomed_size);
    new_way_point = way_points(:,branching_node);
    new_way_point(genes_involved) = 1 - new_way_point(genes_involved);
    way_points = [way_points, new_way_point];
    
    % update terminal nodes
    terminal_nodes = union(setdiff(terminal_nodes, branching_node), size(way_points,2));   
end


% simulate a trajectory-like dataset
num_cells = 3000;
all_data = zeros(num_genes, num_cells);
for i=1:num_cells
    i
    edge_ind = randsample(1:size(directed_edges,1),1);
    mix_alfa = rand;
    cell_expr = way_points(:,directed_edges(edge_ind,1))*(1-mix_alfa) + way_points(:,directed_edges(edge_ind,2))*mix_alfa;
    all_data(:,i) = cell_expr;
end
all_data = all_data + (rand(size(all_data))-0.5)/3;



addpath(genpath('tools\'))
G = digraph(directed_edges(:,1), directed_edges(:,2));
figure(1); subplot(1,3,1)
plot(G)


[coeff, score, latent, tsquared, explained, mu] = pca(all_data');
figure(1); subplot(1,3,2)
scatter(score(:,1),score(:,2),30,'fill');


tSNE_map = fast_tsne(score(:,1:50));
figure(1); subplot(1,3,3)
scatter(tSNE_map(:,1),tSNE_map(:,2),30,'fill');




