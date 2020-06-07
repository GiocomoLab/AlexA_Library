files = dir('/Volumes/Samsung_T5/examples_MEC_speed/*.pdf');
%%
examples = cell(numel(files),1);
fid = fopen('examples_MEC.m','w');
fprintf(fid,'cells = {\n')
for iF=3:numel(files)
    parts = strsplit(files(iF).name,{'_','.'});
    sn = strcat(parts{1},'_',parts{2},'_',parts{3},'_',parts{4},'.mat');
    cluID = str2double(parts{5}(2:end));
    trials = parts{6}(2:end);
    trials = strsplit(trials,{'-',':'});
    trials = str2double(trials{1}):str2double(trials{2});
    examples{iF}={sn,cluID,trials};
    fprintf(fid,'{''%s'',%d,%d:%d},...\n',sn,cluID,trials(1),trials(end));
end
fprintf(fid,'};\n');
fclose(fid)


%%
files = dir('/Volumes/Samsung_T5/examples_V1/*.pdf');
%%
examples = cell(numel(files),1);
for iF=1:numel(files);
    parts = strsplit(files(iF).name,'_');
    sn = strcat(parts{1},'_',parts{2},'_',parts{3},'_',parts{4},'.mat');
    cluID = str2double(parts{5}(2:end));
    trials = parts{8}(3:end);
    trials = strsplit(trials,'-');
    trials = str2double(trials{1}):str2double(trials{2});
    examples{iF}={sn,cluID,trials};
end