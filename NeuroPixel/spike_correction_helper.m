%load matfile data


%OAK = '/oak/stanford/groups/giocomo';
OAK = 'Z:\giocomo';
%NP_DIR = 'F:\NP_DATA';
NP_DIR = fullfile(OAK,'attialex','NP_DATA');
matfiles = dir(fullfile(NP_DIR,'*.mat'));

%dest_path = 'F:/NP_DATA2';
dest_path = fullfile(OAK,'attialex','NP_DATA_corrected');
if ~isfolder(dest_path)
    mkdir(dest_path)
end
uncorrectedList = {};
for iF=1:numel(matfiles)
    if ~isfile(fullfile(dest_path,matfiles(iF).name))
        uncorrectedList{end+1} = matfiles(iF).name;
    end
end
%%

li2={};
for iF=1:numel(uncorrectedList)
    if contains(uncorrectedList{iF},{'gain','contrast','baseline','dark'}) && ~contains(uncorrectedList{iF},{'combined','playbac'})
        li2{end+1}=uncorrectedList{iF};
    end
end 
li2=li2';
%%
% for iF=1:33
%     src = fullfile(NP_DIR,li2{iF});
%     dest = fullfile(OAK,'attialex','NP_DATA_corrected');
%     copyfile(src,dest)
% end