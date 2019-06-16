%% verifying delay and xcorr


pre_id = 10;
post_id = 26;
fr_post = histcounts(sp.st(sp.clu==post_id),0:0.001:max(sp.st));

fr_pre=histcounts(sp.st(sp.clu==pre_id),0:0.001:max(sp.st));
fr_preShift=histcounts(sp.st(sp.clu==pre_id)+0.01,0:0.001:max(sp.st));

[ax,bx]=xcorr(fr_pre,fr_post,100);
[as,bs]=xcorr(fr_pre,fr_preShift,100);
figure
plot(bs/1000,as)
%peak should be at -10, inicating a 10 bin delay in fr_preShift

idx = sp.clu == post_id | sp.clu == pre_id;
[CGR,b]=CCG(sp.st(idx),double(sp.clu(idx))+1,'binSize',[0.001],'duration',[0.2]);


figure
plot(bx/1000,ax)
hold on
plot(b,CGR(:,post_id+1,pre_id+1))
grid on
%this way, CCG output and xcorr output match. i.e. when peak is to the left
%of 0, it means that spiking in pre is followed by spiking in post
 
 %%

postsyn = unique(connected(:,2));
fig = figure;
CGR=mono.ccgR;
for ip=1:length(postsyn)
    n_pairs = nnz(connected(:,2)==postsyn(ip));
    pairs = connected(connected(:,2)==postsyn(ip),1);
    n_rows = 3;
    n_cols = ceil((n_pairs+1)/n_rows);
    
    ha = tight_subplot(n_rows,n_cols,[.01 .03],[.1 .01],[.01 .01]);
    post_idx = postsyn(ip)+1;
    axes(ha(1));
    plot(b,squeeze(CGR(:,post_idx,post_idx)));
    for ic=1:n_pairs
        pre_idx = pairs(ic)+1;

        if mono.Pcausal(pre_idx,post_idx)<0.05
            col = 'r';
        elseif mono.Pcausal(post_idx,pre_idx)<0.05
            col = 'g';
        else
            error('no sig connection')
        end
        axes(ha(ic+1))
        plot(b,squeeze(CGR(:,post_idx,pre_idx)),col);
        grid on
        hold on
        plot(b,squeeze(CGR(:,pre_idx,pre_idx)),'k--')
        xlim([-0.05 0.05])
        %xlabel(sprintf('dist: %.3f',distance(connected(:,1)==pre_idx-1 & connected(:,2)==post_idx-1)))
    end
    pause
    clf;
end
%%
figure
for iP=1:length(connected)
    spike_id=sp.clu==connected(iP,1);
spike_t = sp.st(spike_id);
[~,~,spike_idx] = histcounts(spike_t,post);
subplot(1,2,1)
scatter(posx(spike_idx),trial(spike_idx),2)

subplot(1,2,2)
 spike_id=sp.clu==connected(iP,2);
spike_t = sp.st(spike_id);
[~,~,spike_idx] = histcounts(spike_t,post);

scatter(posx(spike_idx),trial(spike_idx),2)

pause
clf
end

%%

postsyn = unique(connected(:,2));
fig = figure;
for ip=1:length(postsyn)
    n_pairs = nnz(connected(:,2)==postsyn(ip));
    pairs = connected(connected(:,2)==postsyn(ip),1);
    n_rows = 3;
    n_cols = ceil((n_pairs+1)/n_rows);
    
    subplot(n_rows,n_cols,1)
    post_idx = find(good_cells==postsyn(ip));
    plot(DATA.loc,DATA.ACG(post_idx,:));
    for ic=1:n_pairs
        pre_idx = find(good_cells==pairs(ic));
        subplot(n_rows,n_cols,ic+1)
        plot(DATA.loc,DATA.ACG(pre_idx,:));
        grid on
        xlabel(sprintf('dist: %.3f',distance(connected(:,1)==pre_idx-1 & connected(:,2)==post_idx-1)))
    end
    pause
    clf;
end

%%
postsyn = unique(connected(:,2));
fig = figure;
for ip=1:length(postsyn)
    n_pairs = nnz(connected(:,2)==postsyn(ip));
    pairs = connected(connected(:,2)==postsyn(ip),1);

    post_idx = find(good_cells==postsyn(ip));
    plot(DATA.loc,DATA.ACG(post_idx,:));
    hold on
        pre_idx = ismember(good_cells,pairs);
        
        plot(DATA.loc,DATA.ACG(pre_idx,:));
        grid on
   
    pause
    clf;
end
%%
postsyn = unique(connected(:,2));
fig = figure;
for ip=1:length(postsyn)
    n_pairs = nnz(connected(:,2)==postsyn(ip));
    pairs = connected(connected(:,2)==postsyn(ip),1);
    n_rows = 3;
    n_cols = ceil((n_pairs+1)/n_rows);
    
    %subplot(n_rows,n_cols,1)
    ha=tight_subplot(n_rows,n_cols,[.01 .03],[.1 .01],[.01 .01]);
    axes(ha(1));
    post_idx = postsyn(ip)+1;
    spike_id = baseline.sp.clu==postsyn(ip);
    spike_t = baseline.sp.st(spike_id);
[~,~,spike_idx] = histcounts(spike_t,baseline.post);
scatter(baseline.posx(spike_idx),baseline.trial(spike_idx),2)
    for ic=1:n_pairs
        pre_idx = pairs(ic)+1;
        %subplot(n_rows,n_cols,ic+1)
        axes(ha(ic+1));
    spike_id = baseline.sp.clu==pairs(ic);
    spike_t = baseline.sp.st(spike_id);
[~,~,spike_idx] = histcounts(spike_t,baseline.post);
     if mono.Pcausal(pre_idx,post_idx)<0.05
            col = 'r';
        elseif mono.Pcausal(post_idx,pre_idx)<0.05
            col = 'g';
        else
            error('no sig connection')
     end
    scatter(baseline.posx(spike_idx),baseline.trial(spike_idx),2,col)
    end
    pause
    clf;
end
%%
postsyn = unique(connected(:,2));
fig = figure;
for ip=1:length(postsyn)
    n_pairs = nnz(connected(:,2)==postsyn(ip));
    pairs = connected(connected(:,2)==postsyn(ip),1);
    n_rows = 3;
    n_cols = ceil((n_pairs+1)/n_rows);
    
    %subplot(n_rows,n_cols,1)
    ha=tight_subplot(n_rows,n_cols,[.01 .03],[.1 .01],[.01 .01]);
    axes(ha(1));
    post_idx = postsyn(ip)+1;
    imagesc(correlation_All(:,:,post_idx))
    for ic=1:n_pairs
        pre_idx = pairs(ic)+1;
        %subplot(n_rows,n_cols,ic+1)
        axes(ha(ic+1));
    spike_id = baseline.sp.clu==pairs(ic);
    spike_t = baseline.sp.st(spike_id);
[~,~,spike_idx] = histcounts(spike_t,baseline.post);
     if mono.Pcausal(pre_idx,post_idx)<0.05
            col = 'r';
        elseif mono.Pcausal(post_idx,pre_idx)<0.05
            col = 'g';
        else
            error('no sig connection')
     end
    scatter(baseline.posx(spike_idx),baseline.trial(spike_idx),2,col)
    end
    pause
    clf;
end