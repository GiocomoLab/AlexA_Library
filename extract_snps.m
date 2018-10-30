function [snps,aux,trigIDX]=extract_snps(trace,trigs,varargin)
% calculate snps matrix of trig responses per cell
% usage snps=get_trig_response(proj_meta,1,1,trigs,[-50 200]);
% GK 04.11.2016

p = inputParser;
   defaultwin = [-100 100];
   defaultAux=[];
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

% remove any trigger too close to start, finish, or stack transition
toDel=find(sum(abs(bsxfun(@minus,[0 length(trace)]',trigs))<=max(abs(win))));
trigIDX=true(size(trigs));
trigIDX(toDel)=false;
trigs(toDel)=[];

if ~isempty(auxData) && size(AuxData,2) ~= length(trace)
    Error
end

snps=reshape(trace(:,[win(1):win(2)]'*ones(length(trigs),1)'+ones(sum(abs(win))+1,1)*trigs),[size(act,1),sum(abs(win))+1,length(trigs)]);


for iF=1:length(fields)
    if isfield(proj_meta(siteID).rd(masterLayer,tp),fields{iF})
        tmp=proj_meta(siteID).rd(masterLayer,tp).(fields{iF})([win(1):win(2)]'*ones(length(trigs),1)'+ones(sum(abs(win))+1,1)*trigs);
    else
        disp(['no field: ' fields{iF} ', replace with nans'])
        tmp=nan(length(win(1):win(2)),length(trigs));
    end
    aux.(fields{iF})=tmp;
end









