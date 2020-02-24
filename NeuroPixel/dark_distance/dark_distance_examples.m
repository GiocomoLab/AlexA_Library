files = {'Z:\giocomo\attialex\distance_coding3\data/npI1_0414_dark_1.mat',
        'Z:\giocomo\attialex\distance_coding3\data/npH3_0401_dark_3.mat',
        'Z:\giocomo\attialex\distance_coding3\data/npI5_0414_dark_1.mat',
        'Z:\giocomo\attialex\distance_coding3\data/npJ4_0515_dark_1.mat',
        'Z:\giocomo\attialex\distance_coding3\data/npI4_0424_dark_1.mat',
        'Z:\giocomo\attialex\distance_coding3\data/npJ1_0525_dark_1.mat',
        'Z:\giocomo\attialex\distance_coding3\data/npJ5_0504_dark_1.mat'};

    files = files([1 2 3 4])
masterfig = figure;
    
    
for iF=1:numel(files)
    figure
    load(files{iF})
    [~,sid]=sort(mec_depth,'descend');
    v_idx =1./f_vec > 25 &  1./f_vec<600;

    mask = PXX>upper_bound_pxx';
    mask(~v_idx,:)=false;
    fi=medfilt2(double(mask),[3,1]);
    mask=logical(fi);
    subplot(1,4,1)
    imagesc(mask(:,sid)');
    subplot(1,4,2)
    imagesc(ACG(sid,:),[0 0.3])
    title('all')
    [ma,mi]=max(PXX(v_idx,:));
offset = strfind(v_idx,[0 1]);
mi=mi+offset;

keep_PXX=false(size(mi));
i_vec=0;
for ii=1:numel(keep_PXX)
    if all(PXX(mi(ii)+i_vec,ii)>upper_bound_pxx(ii,mi(ii)+i_vec)') && firing_rate(ii)>1 
        keep_PXX(ii)=true;
    end
end
tmp = ACG(keep_PXX,:);
[~,sid]=sort(mec_depth(keep_PXX),'descend');
tmp_PXX = tmp(sid,:);

subplot(1,4,3)
imagesc(tmp_PXX,[0 0.4])
title('PXX Peak based')

keep = false(1,numel(peak_list));
mins = nan(1,numel(peak_list));
for iC=1:numel(peak_list)
       mi = min(ACG(iC,:));
        if ~isempty(peak_list(iC).peak_loc)
            [map,mip]=max(peak_list(iC).peak_val);
            if peak_list(iC).quantile(mip)<map && (map-mi)>.2
                keep(iC)=true;
            end
            mins(iC)=map-mi;
        end
    

end

subplot(1,4,4)
[~,sid2]=sort(mins);
tmp = ACG(keep,:);

td = mec_depth(keep);
[~,sid]=sort(td,'descend');
tmp_XCORR = tmp(sid,:);
imagesc(tmp_XCORR,[0 0.4])
title('xcorr based')

figure(masterfig)
subplot(1,4,iF)

imagesc(tmp_PXX,[0 0.4])
[~,sn,~]=fileparts(files{iF});
title(sn(3:4))
colorbar
end
colormap(summer)
set(gcf,'Renderer','Painters')
