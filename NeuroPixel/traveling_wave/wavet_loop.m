root = 'Y:\giocomo\attialex\NP_DATA';
%Files = dir(fullfile(root,'*_contrast_*.mat'));
Files = dir(fullfile(root,'*gain*.mat'));
MERGED=struct;
cntr=1;
for iF=1:numel(Files)
    clearvars -except root Files iF MERGED cntr
    fn = fullfile(root,Files(iF).name);
    dataset = load(fn);
    try
    quantification_session
    MERGED(cntr).average_triggered = aa_spikes;
    MERGED(cntr).max_depth = maxChan_spikes;
    MERGED(cntr).name = Files(iF).name;
    cntr=cntr+1;
    if mod(iF,7) ~=0
        close(spikefig)
    end
    catch
        warning(strcat('something wrong with ',Files(iF).name))
    end

end

%% scatter of max locations
figure
hold on
for iF=1:numel(MERGED);
    [~,max_loc]=max(MERGED(iF).average_triggered,[],2);
    tmp_d=bins;
    tmp_d=tmp_d-bins(MERGED(iF).max_depth);
    plot(tvec_spikes(max_loc),tmp_d,'b.')
end
%% average spikemats

zz=zeros(2*size(MERGED(iF).average_triggered,1),size(MERGED(iF).average_triggered,2),length(MERGED));
n_bins = numel(bins);
for iF=1:length(MERGED)
    offset = 97-MERGED(iF).max_depth;
    aa_norm = bsxfun(@rdivide,MERGED(iF).average_triggered,sum(MERGED(iF).average_triggered,2));
    zz(offset:offset+n_bins-1,:,iF)=aa_norm;
end
figure
imagesc(flipud(nanmean(zz,3)),[0 0.007])

%set(gca,'YTick',linspace(1,size(spikeMat,1),10),'YTickLabel',round(linspace(max(bins),min(bins),10)))
set(gca,'XTick',linspace(1,numel(tvec_spikes),5),'XTickLabel',linspace(min(tvec_spikes),max(tvec_spikes),5))
set(gca,'YTick',[1 50 99 99+50 194],'YTickLabel',[-max(bins) -max(bins)/2 0 max(bins)/2 max(bins)])
xlabel('Time [s]')
ylabel('Distance from max channel')