ops.towerbins = -41:2:41;
ops.towerbincent = ops.towerbins(1:end-1)*.5 + ops.towerbins(2:end)*.5;
savepath = '/Volumes/Samsung_T5/attialex/python_lr';

%shift_dir = sprintf('Z:/giocomo/attialex/images/xcorrv9/%s_0.80_100',region);
matfiles = dir(fullfile(savepath,'*.mat'));


regions ={'MEC','VISp','RS'};
CorrMat = struct();
SHIFT = struct();
FID = struct();
STAB_BL_Gain = struct();
DEPTH = struct();
ERROR_MLD = struct();
ERROR_CLA =struct();
for iR=1:numel(regions)
    CorrMat.(regions{iR})=[];
    SHIFT.(regions{iR})=[];
    STAB_BL_Gain.(regions{iR}) = [];
    FID.(regions{iR})= [];
    DEPTH.(regions{iR})=[];
    ERROR_MLD.(regions{iR}) = [];
    ERROR_CLA.(regions{iR}) = [];
end


idx_m = triu(true(6),1);
for iF=1:numel(matfiles)
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    
    
%     tmp = data_out.corrMat;
%     tmp_bl = nan(size(tmp,1),1);
%     for iC=1:numel(tmp_bl)
%         tmp_m = squeeze(tmp(iC,1:6,1:6));
%         tmp_bl(iC)=nanmean(tmp_m(idx_m));
%         
%     end
%     tmp_bl_gain = nanmean(nanmean(tmp(:,1:6,7:10),2),3);
%     stab = [tmp_bl';tmp_bl_gain'];
%     sm1=data_out.shiftMat;
%     
%     sm1(ismember(sm1,[-20 20]))=nan;
%     shift1=nanmean(nanmean(sm1(:,1:6,7:10),2),3);
%     shifts = [shift1'];
%     tmp_reg = data_out.region;
%     XC=data_out.corrMat;
    ERROR_T=50;
    pos_bin = discretize(mod(data_out.true_pos,400),ops.xbinedges);
    y = ops.xbincent(pos_bin);
    err_mld = y-ops.xbincent(data_out.predicted_bin);
    err_cla = y-ops.xbincent(data_out.predicted_theta+1);
    correction_mld = abs(err_mld)>400/2;
    err_mld(correction_mld) = err_mld(correction_mld)-400*sign(err_mld(correction_mld));
    err_mld(abs(err_mld)>ERROR_T)=nan;
    err_mld=abs(err_mld);
    
    correction_cla = abs(err_cla)>400/2;
    err_cla(correction_cla) = err_cla(correction_cla)-400*sign(err_cla(correction_cla));
    err_cla(abs(err_cla)>ERROR_T)=nan;
    err_cla = abs(err_cla);
    
    avg_cla = zeros(1,numel(ops.nBins));
    avg_mld = avg_cla;
%     if ~isfield(data_out,'trials')
%         tmp = unique(data_out.trial);
%         data_out.decoder.trials = tmp;
%     end
%    trial_idx = ismember(data_out.decoder.trial,data_out.decoder.trials(7:10));
    trial_idx = true(size(pos_bin));
    for ii=1:200
        idx = pos_bin==ii & trial_idx;
        avg_cla(ii)=nanmean(err_cla(idx));
        avg_mld(ii)=nanmean(err_mld(idx));
    end
    
    
    
    counts = nan(size(regions));
    for iR=1:3
        counts(iR)=nnz(startsWith(data_out.region,regions{iR}));
    end
    [a,ii]=max(counts);
    winner = regions{ii};
    if startsWith(winner,'RS')
        winner = 'RS';
    end
    reg_idx = startsWith(data_out.region,winner);
    if nnz(reg_idx)>=1
        
            %CorrMat.(winner) = cat(1,CorrMat.(winner),(nanmean(XC(reg_idx,:,:),1)));
            %SHIFT.(winner) = cat(2,SHIFT.(winner),shifts(:,reg_idx));
            
            FID.(winner) = cat(2,FID.(winner),iF*ones(1,nnz(reg_idx)));
            %DEPTH.(regions{iR}) = cat(2,DEPTH.(regions{iR}),data_out.depth(reg_idx))*mult;
%            STAB_BL_Gain.(winner) = cat(2,STAB_BL_Gain.(winner),stab(:,reg_idx));
            ERROR_MLD.(winner) = cat(1,ERROR_MLD.(winner),avg_mld);
            ERROR_CLA.(winner) = cat(1,ERROR_CLA.(winner),avg_cla);
    end

    end
    
    
    
    
    


%% get cluster1 sessions
ALLC = [];
ALLSHIFTS = [];
ALLR = [];
E_MLD=[];
E_CLA=[];
TE_MLD=[];
TE_CLA = [];
for iR=1:3
    
            %ALLC=cat(1,ALLC,CorrMat.(regions{iR}));
            %ALLSHIFTS = cat(2,ALLSHIFTS,(SHIFT.(regions{iR})));
            ALLR = cat(1,ALLR,iR*ones(size(ERROR_MLD.(regions{iR}),1),1));
            tmp_mld = ERROR_MLD.(regions{iR});
            tmp_cla = ERROR_CLA.(regions{iR});
            E_MLD = cat(1,E_MLD,tmp_mld);
            E_CLA = cat(1,E_CLA,tmp_cla);
            t1 = nanmean(cat(3,tmp_mld(:,20:60),tmp_mld(:,60:100),tmp_mld(:,100:140),tmp_mld(:,140:180)),3);
            %t1 = nanmean(cat(3,tmp_mld(:,20:60),tmp_mld(:,60:100)),3);
            TE_MLD = cat(1,TE_MLD,t1);
            t2=nanmean(cat(3,tmp_cla(:,20:60),tmp_cla(:,60:100),tmp_cla(:,100:140),tmp_cla(:,140:180)),3);
            %t2=nanmean(cat(3,tmp_cla(:,20:60),tmp_cla(:,60:100)),3);
            TE_CLA = cat(1,TE_CLA,t2);
            
            
        
    
end

tf = triu(true(16),1);
 
% [coeff,score,~,~,expl] = pca(ALLC(:,tf));
% 
% % K means clustering
% rng(2); % for reproducibility
% rep_cluster = kmeans(score(:,1:3),3);
% 
% figure('Position',[680         817        1080         281])
% scatter(score(:,1),score(:,2),15,rep_cluster)
% cluster2look = 1;
% figure
% for ii=1:3
%     IDX = ALLR == ii;
% 
%     subplot(1,3,ii)
%     hold on
%     imagesc(1:16,1:15,squeeze(nanmean(ALLC(IDX,:,:),1)))
% end

figure('Position',[680         511        1019         587])
cmap=cbrewer('qual','Set2',3,'pchip');

for ii=1:3
        IDX = ALLR == ii;
    subplot(2,3,ii)
    nn = nanmean(E_MLD(IDX,:));
    ee = nanstd(E_MLD(IDX,:))/sqrt(nnz(IDX));

    %plot(ops.xbincent,tmp)
    boundedline(ops.xbincent,nn,ee,'alpha','cmap',[0 0 1])
    
    nn = nanmean(E_CLA(IDX,:));
    ee = nanstd(E_CLA(IDX,:))/sqrt(nnz(IDX));
    boundedline(ops.xbincent,nn,ee,'alpha','cmap',[1 0 0])
        xlabel('Distance in tunnel')
    title(regions{ii})
    legend('MLD Decoder','Classifier')
    


    subplot(2,3,ii+3)
    nn = nanmean(TE_MLD(IDX,:));
    ee = nanstd(TE_MLD(IDX,:))/sqrt(nnz(IDX));

    %plot(ops.xbincent,tmp)
    boundedline(ops.towerbincent,nn,ee,'alpha','cmap',[0 0 1])
    
        nn = nanmean(TE_CLA(IDX,:));
    ee = nanstd(TE_MLD(IDX,:))/sqrt(nnz(IDX));

    %plot(ops.xbincent,tmp)
    boundedline(ops.towerbincent,nn,ee,'alpha','cmap',[1 0 0])
    xlabel('Distance from towers')
    title(regions{ii})
        legend('MLD Decoder','Classifier')

end