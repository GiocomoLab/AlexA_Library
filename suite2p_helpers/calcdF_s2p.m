function [ dF ] = calcdF_s2p(input_dir)

F = readNPY(fullfile(input_dir,'F.npy'));
Fneu = readNPY(fullfile(input_dir,'Fneu.npy'));
iscell = readNPY(fullfile(input_dir,'iscell.npy'));
cellIDX=iscell(:,1)==1;
%%
dF=F(cellIDX,:)-0*Fneu(cellIDX,:);
%%
% win=gausswin(20)';
% win=win/sum(win);
% Flow = conv2(dF,win,'same');
% 
% Flow = imerode(Flow,ones(1,60));
% Flow = imdilate(Flow,ones(1,60));
% 
% dF=dF-Flow;

%%
nRois=size(dF,1);
    for gnd=1:nRois
        cur_act=psmooth(dF(gnd,:),30,15);
        dF(gnd,:)=cur_act/median(cur_act);
    end


end

