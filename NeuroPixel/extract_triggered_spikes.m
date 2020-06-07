function [spike_mat,win,aux_mat,mean_response] = extract_triggered_spikes(sp_struct,time_idx, varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
p = inputParser;
defaultwin = [-.5 3];
defaultaux = [];
defaultAuxWin = [-25 150];
addParameter(p,'win',defaultwin);
addParameter(p,'aux',defaultaux); % assumes that first row is time vec
addParameter(p,'aux_win',defaultAuxWin)

parse(p,varargin{:});

win=p.Results.win;
aux=p.Results.aux;
aux_win=p.Results.aux_win;

spike_times=sp_struct.st;

time_idx(time_idx+win(1)<=0)=[];
time_idx(time_idx+win(2)>max(spike_times))=[];
n_times = length(time_idx);
cluster_ID=sp_struct.clu;
n_units = max(cluster_ID)+1;
tmp = max(sp_struct.cids)+1;
n_units=max(n_units,tmp);

n_aux = size(aux,1)-1;

time_vec=win(1):0.001:win(2);
spike_mat=zeros(n_units,n_times,length(time_vec));

for iT=1:n_times
    start = time_idx(iT)+win(1);
    stop = time_idx(iT)+win(2);
    spike_idx = (spike_times>start & spike_times<stop);
    cluIDs = cluster_ID(spike_idx)+1;
    times = round((spike_times(spike_idx)-start)*1000)+1;
    IND = sub2ind(size(spike_mat),cluIDs,iT*ones(size(times)),times);
    spike_mat(IND)=1;
    
    
end
spike_mat = spike_mat(sp_struct.cids+1,:,:);

if ~isempty(aux)
    idx_win=aux_win(1):aux_win(2);
    aux_mat = zeros(n_aux,n_times,length(idx_win));
    aux_time=aux(1,:);
    for iT=1:n_times
        %find idx closest to time
        [~,idx]=min(abs(aux_time-time_idx(iT)));
        aux_IDX=idx+idx_win;
        aux_mat(:,iT,:)=aux(2:end,aux_IDX);
    end
end

if nargout == 4
    kernel=reshape(gausswin(401),1,1,[]);
    binsize=0.001;
    kernel=kernel/sum(kernel(:))/binsize;
    rate=convn(spike_mat,kernel,'same');
    mean_response=squeeze(mean(rate,2));
end