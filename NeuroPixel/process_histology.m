function histology = process_histology()
if ispc()
    histo_table = readtable('Z:\giocomo\attialex\NP_DATA\Histology_Quantification.xlsx');
elseif ismac()
    histo_table = readtable('/Users/attialex/code/AlexA_Library/NeuroPixel/histology/Histology_Quantification.xlsx');
else
    histo_table = readtable('/oak/stanford/groups/giocomo/attialex/NP_DATA/Histology_Quantification.xlsx');
end

time_dat = cell(numel(histo_table.Date),1);
for iD=1:numel(histo_table.Date);
    time_dat{iD}=datestr(histo_table.Date(iD),'mmdd');
end

year = datestr(histo_table.Date,'yy');
histology.date = time_dat;
histology.year = year;
histology.animal = histo_table.Mouse;
n_entries = size(histo_table,1);
origin = zeros(n_entries,3);
unit_vec = zeros(n_entries,3);
for iT=1:n_entries
mec_entry = [histo_table.X1(iT),histo_table.Y1(iT),histo_table.Z1(iT)];
probe_term = [histo_table.X2(iT),histo_table.Y2(iT),histo_table.Z2(iT)];

vec1 = mec_entry-probe_term;
    vec1(3) = -vec1(3); % flip Z coord
    vec1 = vec1./repmat(sqrt(sum(vec1.^2,2)),1,3); % normalize
    vec2 = probe_term; 
    vec2(3) = -vec2(3); % flip Z coord
    unit_vec(iT,:) = vec1;
    origin(iT,:) = vec2;
    
end
histology.unit_vector = unit_vec;
histology.origin = origin;
for iF=1:length(histo_table.Properties.VariableNames)
    var= histo_table.Properties.VariableNames{iF};
    histology.(var) = histo_table.(var);
end

end