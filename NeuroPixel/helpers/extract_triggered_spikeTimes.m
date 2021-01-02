function [spike_times_struct,win,aux_mat,time_idx_used,count_vec] = extract_triggered_spikeTimes(sp_struct,time_idx, varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
p = inputParser;
defaultwin = [-.5 3];
defaultaux = [];
defaultAuxWin = [-25 150];
addParameter(p,'win',defaultwin);
addParameter(p,'aux',defaultaux); % assumes that first row is time vec
addParameter(p,'aux_win',defaultAuxWin)
addParameter(p,'cluIDs',[])
addParameter(p,'opt',[])

parse(p,varargin{:});

win=p.Results.win;
aux=p.Results.aux;
aux_win=p.Results.aux_win;

if isempty(p.Results.cluIDs)
    %only extract good cells
    good_cells = sp_struct.cids(sp_struct.cgs==2);
    take_idx = ismember(sp_struct.clu,good_cells);
else
    take_idx = ismember(sp_struct.clu,p.Results.cluIDs);
    good_cells = p.Results.cluIDs;
end

if isempty(p.Results.opt)
    opt = load_mismatch_opt;
end


spike_times=sp_struct.st(take_idx);
take_idx_time = true(size(time_idx));
take_idx_time(time_idx+win(1)<=0)=false;
take_idx_time(time_idx+win(2)>max(spike_times))=false;

time_idx(time_idx+win(1)<=0)=[];
time_idx(time_idx+win(2)>max(spike_times))=[];
time_idx_used = take_idx_time;
n_times = length(time_idx);
cluster_ID=sp_struct.clu(take_idx);


n_aux = size(aux,1)-1;

spike_times_struct = cell(n_times,1);

for iT=1:n_times
    start = time_idx(iT)+win(1);
    stop = time_idx(iT)+win(2);
    spike_idx = (spike_times>start & spike_times<stop);
    cluIDs = cluster_ID(spike_idx);
    times = (spike_times(spike_idx)-start+win(1));
    spike_times_struct{iT}=cat(2,times,double(cluIDs),ones(size(times))*iT);
end

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
    
else
    aux_mat = [];
end

trial_vec =cat(1,spike_times_struct{:});
count_vec = zeros(numel(good_cells),numel(opt.time_bins)-1);
for iC=1:numel(good_cells)
    idx = trial_vec(:,2)==good_cells(iC);
    [spike_count]=histcounts(trial_vec(idx,1),opt.time_bins);
    count_vec(iC,:)=spike_count;
    
end
n_trigs_included = numel(unique(trial_vec(:,3)));
count_vec = count_vec/n_trigs_included/opt.TimeBin;

