if ispc()
    data_dir = 'Z:\giocomo\attialex\NP_DATA\';
    image_dir = 'Z:\giocomo\attialex\images\dark_acg\';
else
    data_dir = '/oak/stanford/groups/giocomo/attialex/NP_DATA';
    image_dir = '/oak/stanford/groups/giocomo/attialex/images';
end
h=figure('Position',[  446         281        1375         697]);
files = dir(fullfile(data_dir,'*dark*.mat'));
params = readtable('UniversalParams.xlsx');
for iF=11:length(files)
    clear DATA
    load(fullfile(data_dir,files(iF).name));
    
    idx = (DATA.loc*2)>25 & (DATA.loc*2)<2000;
    start = find(idx,1);
    start=start-1;
    good_cells = sp.cids(sp.cgs==2);
    
    spacing = DATA.loc;
    av_speed = total_distance(end)/post(end);
    for ii=1:length(good_cells)
        ACG=DATA.ACG(ii,:);
        subplot(2,3,1)
        plot(DATA.loc*2,ACG)
        hold on
        plot(spacing*2,DATA.ACG_RANDOM(ii,:),'r')
        plot(spacing*2,DATA.ACG_QUANTILES{ii}','r--')
        xlim([-2000 2000])
        xlabel('cm')
        legend('ACG','Shuffle','.9 Quantile')
        [a,b]=sort(ACG(idx),'descend');
        diff_left = a-ACG(b+start-1);
        diff_right = a-ACG(b+start+1);
        
        pot_idx = find(diff_left>0 & diff_right > 0,1);
        has_min = nnz(ACG<DATA.ACG_QUANTILES{ii}(2,:))>0;
        if ~isempty(pot_idx) && has_min 
            pot_idx = b(pot_idx(1))+start;
            above_noise = ACG(pot_idx)>DATA.ACG_QUANTILES{ii}(1,pot_idx);
            
            if above_noise
                plot(spacing(pot_idx)*2,ACG(pot_idx),'ro')
            end
            
        end
        title('Spatial ACG')
        subplot(2,3,4)
        spike_id = sp.clu == good_cells(ii);
        if ~isempty(pot_idx) && nnz(spike_id)>100
            
            [CGR,b]=CCG(sp.st(spike_id),double(sp.clu(spike_id)+1),'binSize',[0.1],'duration',[80]);
            plot(b,CGR(:,good_cells(ii)+1,good_cells(ii)+1))
            if ~isempty(pot_idx) && has_min && above_noise
                x=spacing(pot_idx)*2/av_speed;
                y=interp1(b,CGR(:,good_cells(ii)+1,good_cells(ii))+1,x);
                hold on
                plot(x,y,'ro')
            end
        end
        title('temporal ACG')
        xlabel('time [s]')
        subplot(2,3,[2 5])
        spike_t = sp.st(spike_id);
        [~,~,spike_idx] = histcounts(spike_t,post);
        
        scatter(posx(spike_idx),trial(spike_idx),2)
        axis tight
        
        if ~isempty(pot_idx) && has_min && above_noise && nnz(spike_id)>100
            
            %plot(spacing(pot_idx)*2,ACG(pot_idx),'ro')
            subplot(2,3,[3 6])
            cell_spacing = spacing(pot_idx)*2;
            trial_tmp=floor(total_distance/cell_spacing);
            pos_tmp = mod(total_distance,cell_spacing);
            
            scatter(pos_tmp(spike_idx),trial_tmp(spike_idx),2)
            axis tight
            
            
            %pause
            saveas(h,fullfile(image_dir,sprintf([files(iF).name '_%d.png'],good_cells(ii))),'png');
        end
        clf
    end
end