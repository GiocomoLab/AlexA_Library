function [snps,aux,trigIDX]=extract_snps(trace,trigs,varargin)
% calculate snps matrix of trig responses per cell
% usage snps=extract_snps(cell_activity,trigs,'win',[-50 200]);
% cell_activity: neurons x frames, trigs: frame idx of onset, win: extracts frames from win(1) to win(2) around trigger
% 'aux': a struct containing aux data to extract. E.g. adata.running, adata.stimID (1 x frames)

p = inputParser;
   defaultwin = [-100 100];
   defaultAux=struct();
   %defaultLayer = 1:4;
   %defaultFields = {'velM_smoothed','velP_smoothed'};

   addParameter(p,'win',defaultwin);
   addParameter(p,'aux',defaultAux);
   %addParameter(p,'layer',defaultLayer);
   %addParameter(p,'fields',defaultFields);
  
   
   parse(p,varargin{:});

trigs=trigs(:)';

win=p.Results.win;
auxData=p.Results.aux;
fields = fieldnames(auxData);

% remove any trigger too close to start, finish, or stack transition
toDel=find(sum(abs(bsxfun(@minus,[0 length(trace)]',trigs))<=max(abs(win))));
trigIDX=true(size(trigs));
trigIDX(toDel)=false;
trigs(toDel)=[];


snps=reshape(trace(:,[win(1):win(2)]'*ones(length(trigs),1)'+ones(sum(abs(win))+1,1)*trigs),[size(trace,1),sum(abs(win))+1,length(trigs)]);


for iF=1:length(fields)
    if isfield(auxData,fields{iF})
        tmp=auxData.(fields{iF})([win(1):win(2)]'*ones(length(trigs),1)'+ones(sum(abs(win))+1,1)*trigs);
    else
        disp(['no field: ' fields{iF} ', replace with nans'])
        tmp=nan(length(win(1):win(2)),length(trigs));
    end
    aux.(fields{iF})=tmp;
end









