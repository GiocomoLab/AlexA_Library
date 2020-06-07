
savepath = '/Volumes/Samsung_T5/tbtxcorr_with_shortbaseline2';

%shift_dir = sprintf('Z:/giocomo/attialex/images/xcorrv9/%s_0.80_100',region);
matfiles = dir(fullfile(savepath,'*.mat'));


regions ={'MEC','VISp','RS'};
CorrMat = struct();
SHIFT = struct();
FID = struct();
STAB_BL_Gain = struct();
DEPTH = struct();
for iR=1:numel(regions)
    CorrMat.(regions{iR})=[];
    SHIFT.(regions{iR})=[];
    STAB_BL_Gain.(regions{iR}) = [];
    FID.(regions{iR})= [];
    DEPTH.(regions{iR})=[];
end


idx_m = triu(true(6),1);
for iF=1:numel(matfiles)
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    
    
    tmp = data_out.corrMat;
    tmp_bl = nan(size(tmp,1),1);
    for iC=1:numel(tmp_bl)
        tmp_m = squeeze(tmp(iC,1:6,1:6));
        tmp_bl(iC)=nanmean(tmp_m(idx_m));
        
    end
    tmp_bl_gain = nanmean(nanmean(tmp(:,1:6,7:10),2),3);
    stab = [tmp_bl';tmp_bl_gain'];
    sm1=data_out.shiftMat;
    sm2 =data_out.shiftMatShifted;
    sm1(ismember(sm1,[-20 20]))=nan;
    sm2(ismember(sm2,[-20 20])) = nan;
    shift1=nanmean(nanmean(sm1(:,1:6,7:10),2),3);
    shift2=nanmean(nanmean(sm2(:,1:6,7:10),2),3);
    shifts = [shift1';shift2'];
    tmp_reg = data_out.region;
    XC=data_out.corrMat;
    
    for iR=1:numel(regions)
        reg_idx = startsWith(tmp_reg,regions{iR});
        if ismember(regions{iR},{'MEC','ECT'})
            mult = -1;
        else
            mult = 1;
        end
        if any(reg_idx)
            CorrMat.(regions{iR}) = cat(1,CorrMat.(regions{iR}),XC(reg_idx,:,:));
            SHIFT.(regions{iR}) = cat(2,SHIFT.(regions{iR}),shifts(:,reg_idx));
            
            FID.(regions{iR}) = cat(2,FID.(regions{iR}),iF*ones(1,nnz(reg_idx)));
            DEPTH.(regions{iR}) = cat(2,DEPTH.(regions{iR}),data_out.depth(reg_idx))*mult;
            STAB_BL_Gain.(regions{iR}) = cat(2,STAB_BL_Gain.(regions{iR}),stab(:,reg_idx));
        end
    end
    
    
end
%% plot across cells
figure
X=[];
G=[];
TMP = {};
for iR=1:3
    stab = STAB_BL_Gain.(regions{iR})(2,:);
    idx = stab>.7;
    shifts = SHIFT.(regions{iR});
    shifts(ismember(shifts,[-20 20]))=nan;
    subplot(2,3,iR)
    histogram(shifts(1,idx),[-20:20]);
    hold on
    histogram(shifts(2,idx),[-20:20]);
    subplot(2,3,iR+3)
    scatter(shifts(1,idx),shifts(2,idx))
    hold on
    plot([-30 30],[-30 30],'k')
    grid on
    axis image
    xlim([-21 21])
    ylim([-21 21])
    X=cat(2,X,diff(shifts(:,idx)));
    G=cat(2,G,iR*ones(1,nnz(idx)));
    TMP.(regions{iR})=diff(shifts(:,idx));
end
%% get cluster1 sessions
ALLC = [];
ALLSHIFTS = [];
ALLR = [];
for iR=1:3
    sites = FID.(regions{iR});
    usites=unique(sites);
    for iS=1:numel(usites)
        sid = sites==usites(iS);
        stab_idx = STAB_BL_Gain.(regions{iR})(1,:)>.5;
        if nnz(sid & stab_idx)>=5
            idx = sid & stab_idx;
            ALLC=cat(1,ALLC,nanmean(CorrMat.(regions{iR})(idx,:,:)));
            ALLSHIFTS = cat(2,ALLSHIFTS,nanmean(SHIFT.(regions{iR})(:,idx),2));
            ALLR = cat(1,ALLR,iR);
        end
    end
end

tf = triu(true(16),1);
 
[coeff,score,~,~,expl] = pca(ALLC(:,tf));

% K means clustering
rng(2); % for reproducibility
rep_cluster = kmeans(score(:,1:3),3);

figure('Position',[680         817        1080         281])
% subplot(1,3,1)
% scatter(score(:,1),score(:,2),15,rep_cluster)
subplot(1,3,1)
hold on
mt = {'o','x','v'};
for iR=1:3
 scatter(ALLSHIFTS(1,ALLR==iR & rep_cluster==1),ALLSHIFTS(2,ALLR==iR & rep_cluster==1),25,[0 0 0],mt{iR})
end
plot([-15 5],[-15 5],'--','Color',[0 0 0])
legend(regions)
legend('Location','SouthEast')

dd={};
X=[];
G=[];
for iR=1:3
    tmp = ALLSHIFTS(1,ALLR==iR & rep_cluster==1)-ALLSHIFTS(2,ALLR==iR & rep_cluster==1);
 dd{iR}=tmp;
 X=cat(2,X,tmp);
 G=cat(2,G,iR*ones(size(tmp)));
end
subplot(1,3,2)
plotSpread(dd)
title('delta shift')
xticklabels(regions)
p=anova1(X,G,'off');

xlabel(sprintf('anova: %.3e',p))
for ii=1:3
    bar(ii,median(X(G==ii)),'EdgeColor','b','FaceColor','none')

    text(ii-.5,1,sprintf('%.3e',signrank(X(G==ii))))
end
text(1,-3,sprintf('%.3e',ranksum(X(G==1),X(G==2))));
text(2,-3,sprintf('%.3e',ranksum(X(G==1),X(G==2))));
box off
subplot(1,3,3)
S1=[];
S2=[];
G=[];
for iR=1:3
    tmp1 = ALLSHIFTS(1,ALLR==iR & rep_cluster==1);
    tmp2 = ALLSHIFTS(2,ALLR==iR & rep_cluster==1);
    dd{iR}=tmp2;
    S1 = cat(2,S1,tmp1);
    S2 = cat(2,S2,tmp2);
    G = cat(2,G,iR*ones(size(tmp1)));
end
plotSpread(dd)
p=anova1(S2,G,'off');
xticklabels(regions)
xlabel(sprintf('anova: %.3e',p))
box off
for ii=1:3
    bar(ii,median(S2(G==ii)),'EdgeColor','b','FaceColor','none')
    text(ii-.5,2,sprintf('%.3e',signrank(S2(G==ii))))
end
title('Residual Shift')

saveas(gcf,'/Volumes/Samsung_T5/attialex/images/gain_shift_summary.pdf')