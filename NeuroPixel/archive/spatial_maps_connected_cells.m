%% load the data
root = 'Z:\giocomo\attialex\NP_DATA\';
baseline = load(fullfile(root,'npI4_0421_baseline_1'));
dark = load(fullfile(root,'npI4_0421_dark_1'));

%% compute spatial similarity
[sim,sim_s]=compute_spatial_similarity(baseline);
%%
connected = dark.connected;
postsyn = unique(connected(:,2));
fig = figure;
good_cells=baseline.sp.cids(baseline.sp.cgs==2);
for ip=1:length(postsyn)
    n_pairs = nnz(connected(:,2)==postsyn(ip));
    pairs = connected(connected(:,2)==postsyn(ip),1);
    n_rows = 3;
    n_cols = ceil((n_pairs+1)/n_rows);
    
    %subplot(n_rows,n_cols,1)
    ha=tight_subplot(n_rows,n_cols,[.01 .03],[.1 .01],[.01 .01]);
    axes(ha(1));
    post_idx = find(good_cells == postsyn(ip));
    imagesc(sim_s(:,:,post_idx))
    for ic=1:n_pairs
        pre_idx = find(good_cells == pairs(ic));
        %subplot(n_rows,n_cols,ic+1)
        axes(ha(ic+1));
    imagesc(sim_s(:,:,pre_idx));
    end
    pause
    clf;
end