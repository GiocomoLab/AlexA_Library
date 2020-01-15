files = dir('Z:\giocomo\attialex\decodeBayesv3\*.mat');

errors = containers.Map;
errors_deg = containers.Map;
posteriors = containers.Map;
plotfigs = true;
if plotfigs
    figure
end


for iF=1:numel(files)%[6 25]%1:numel(files)
    data = load(fullfile(files(iF).folder,files(iF).name));
    set(gcf,'Name',files(iF).name)
    nR=numel(data.data_out.errors);
    for iR=1:nR
        if plotfigs
            subplot(3,nR,iR)
        end
        if data.data_out.n_cells{iR}
            bin_centers = .5*data.data_out.edges_all{iR}(1:end-1)+.5*data.data_out.edges_all{iR}(2:end);
            if plotfigs
                plot(data.data_out.location_real{iR})
                hold on
                plot(data.data_out.location_decoded{iR});
                title(data.data_out.regions{iR})
                subplot(3,nR,nR+iR)
                imagesc(log(data.data_out.posteriors{iR}))
            end
            Bins = 1:numel(data.data_out.edges_all{iR})-1;
            tr = data.data_out.trial{iR};
            error_tmp = nan(4,numel(Bins));
            error_tmp_deg = error_tmp;
            posterior_tmp = zeros(numel(Bins),numel(Bins),5);
            error = data.data_out.errors{iR};
            error(error>200)=error(error>200)-400;
            error(error<-200) = error(error<-200)+400;
            error = sqrt(error.^2);
            if plotfigs
                subplot(3,nR,2*nR+iR)
                plot(error)
            end
            d_loc = discretize(data.data_out.location_real{iR},bin_centers);
            
            traj_r = data.data_out.location_real{iR};
            traj_d = data.data_out.location_decoded{iR};
            traj_d = traj_d/400*2*pi;
            traj_r = traj_r/400*2*pi;
            radius = 1;
            
            x_r=radius*cos(traj_r);
            y_r=radius*sin(traj_r);
            x_d=radius*cos(traj_d);
            y_d=radius*sin(traj_d);
            v_r=[x_r;y_r;zeros(size(x_r))];
            v_d=[x_d;y_d;zeros(size(x_r))];
            for ii=1:numel(traj_r)
                c=cross(v_r(:,ii),v_d(:,ii));
                d=dot(v_r(:,ii),v_d(:,ii));
                mult = sign(c(3));
                th(ii)=mult*atan2d(norm(c),d);
            end
            
            speed = data.data_out.speed_all{iR};
            thresh =1;
            for iT=1:5
                %disp(data.data_out.gains{iR}(iT))
                for iB=1:numel(Bins)
                    idx = tr==(20+iT) & d_loc==iB & speed>thresh;
                    tmp=nanmean(error(idx));
                    error_tmp(iT,iB)=tmp;
                    tmp_deg = nanmean(th(idx));
                    error_tmp_deg(iT,iB)=tmp_deg;
                    
                    if nnz(idx)==1
                        tmp = data.data_out.posteriors{iR}(:,idx);
                    elseif nnz(idx)>1
                        tmp = nanmean(data.data_out.posteriors{iR}(:,idx),2);
                    else
                        tmp = nan(numel(Bins),1);
                    end
                    posterior_tmp(:,iB,iT)=tmp;
                    
                end
            end
            key = data.data_out.regions{iR};
            if errors.isKey(key)
                tmp = errors(key);
                tmp = cat(1,tmp,error_tmp);
                tmp_deg = errors_deg(key);
                tmp_deg = cat(1,tmp_deg,error_tmp_deg);
                errors(key)=tmp;
                errors_deg(key)=tmp_deg;
                tmp_post = posteriors(key);
                tmp_post = cat(3,tmp_post,posterior_tmp);
                posteriors(key)=tmp_post;
            else
                errors(key)=error_tmp;
                errors_deg(key)=error_tmp_deg;
                posteriors(key)=posterior_tmp;
            end
            
        end
    end
    if plotfigs
    pause
    clf
    end
end

%%
kk=errors.keys();
kk={'MEC'};
for iK=1:numel(kk)
    tmp = errors(kk{iK});
    figure
    plot(bin_centers,nanmean(tmp(1:5:end,:)),'.-')
    hold on
    plot(bin_centers,nanmean(tmp(4:5:end,:)),'.-')
    %plot(bin_centers,nanmean(tmp(3:4:end,:)),'.-')
    %plot(bin_centers,nanmean(tmp(5:5:end,:)),'.-')
    
    
    title(kk{iK})
end
%%
kk=errors.keys();

for iK=1:numel(kk)
    tmp = errors(kk{iK});
    figure
    plot(bin_centers,nanmean(tmp(1:4:end,:)),'.-')
    hold on
    plot(bin_centers,nanmean(tmp(2:4:end,:)),'.-')
    %plot(bin_centers,nanmean(tmp(3:4:end,:)),'.-')
    %plot(bin_centers,nanmean(tmp(4:4:end,:)),'.-')
    
    
    title(kk{iK})
end
%%
meanE = errors('MEC');
meanE = nanmean(meanE,2);
figure
hold on

for iT=1:4
    histogram(meanE(iT:5:end),'DisplayStyle','stairs')
end

%%
meanE = errors('MEC');
meanE = nanmean(meanE,2);
xvec = [1:5];
xvec = repmat(xvec,1,59);
figure
beeswarm(xvec',meanE)
%%
for iT=1:4
    
tmp = posteriors('MEC');
figure
imagesc(log(nanmean(tmp(:,:,1:5:end),3)))
hold on

plot([1,200],[1,200],'r-')
end