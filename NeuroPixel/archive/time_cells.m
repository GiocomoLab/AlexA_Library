%load data spike data

%load running data

%find indices of stop running that are followed by at least 3s below
%threshold 
run_bi = running<run_t;
runwin=[ones(1,4*sr) zeros(1,3*sr)];
stop_idx=strfind(run_bi,runwin)+4*sr;

spikes=extract_triggered_spikes(data,stop_idx,'win',[-100 100]);

%convert to Firing Rate

fr=...
    
%average


%sort by peak


%plot