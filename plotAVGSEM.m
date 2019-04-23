function plotAVGSEM(avg,ax,varargin)
%plots mean and SEM across neurons. Time x Neurons
%   Detailed explanation goes here

params=struct();
params.winIDX=-100:100;
params.masterTime=params.winIDX/15;
params.xLim=[-1 4];

avg=double(avg);


p = inputParser;
   defaultLinetype = '-';
   defaultMeanSubtract = true;
   defaultGroup=[];
   defaultBaseline=[93:100];
  
   defaultPars = params;
   defaultColor = [0 0 1];
   
   addParameter(p,'parameters',defaultPars);
   addParameter(p,'lt',defaultLinetype)
   addParameter(p,'col',defaultColor);
   addParameter(p,'ms',defaultMeanSubtract);
   addParameter(p,'grouping',defaultGroup);
   addParameter(p,'baseline',defaultBaseline);

   parse(p,varargin{:});

    


msSTD=nanstd(avg,[],2)/sqrt(size(avg,2));
params=p.Results.parameters;
params.baseline=p.Results.baseline;

if isempty(p.Results.grouping)
    msAVG=nanmean(avg,2);
    
else
    [a]=unique(p.Results.grouping);
    AVG=zeros(size(avg,1),length(a));
    for ii=1:length(a)
        idx=p.Results.grouping==a(ii);
        AVG(:,ii)=nanmean(avg(:,idx),2);
    end
    msAVG=nanmean(AVG,2);
    msSTD=nanstd(AVG,[],2)/sqrt(length(a));
end


if p.Results.ms
msAVG=msAVG-nanmean(msAVG(params.baseline,:));
end

    


axes(ax)
hold on

 
    
boundedline(params.masterTime,msAVG,msSTD,'cmap',p.Results.col,p.Results.lt,'alpha')

xlim(params.xLim)
end

