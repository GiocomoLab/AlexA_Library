files = {'Z:\giocomo\attialex\distance_coding2\data/npI1_0414_dark_1.mat',
        'Z:\giocomo\attialex\distance_coding2\data/npH3_0401_dark_3.mat',
        'Z:\giocomo\attialex\distance_coding2\data/npI5_0414_dark_1.mat',
        'Z:\giocomo\attialex\distance_coding2\data/npJ4_0515_dark_1.mat',
        'Z:\giocomo\attialex\distance_coding2\data/npJ1_0525_dark_1.mat',
        'Z:\giocomo\attialex\distance_coding2\data/npJ5_0504_dark_1.mat'};

masterfig = figure;
    
    
for iF=1%:numel(files)
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
    
    [ma,mi]=max(PXX(v_idx,:));
offset = strfind(v_idx,[0 1]);
mi=mi+offset;

keep=false(size(mi));
i_vec=-1:1;
for ii=1:numel(keep)
    if all(PXX(mi(ii)+i_vec,ii)>upper_bound_pxx(ii,mi(ii)+i_vec)')
        keep(ii)=true;
    end
end
subplot(1,4,3)
imagesc(ACG(sid(keep),:),[0 0.4])
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
imagesc(tmp(sid,:),[0 0.4])

figure(masterfig)
subplot(1,4,iF)
imagesc(tmp(sid,:),[0 0.4])
[~,sn,~]=fileparts(files{iF});
title(sn(3:4))
colorbar
end
colormap(summer)
set(gcf,'Renderer','Painters')
