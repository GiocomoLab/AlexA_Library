 


function medianTrace = applyCARtoDat_malcolm(filename, nChansTotal, outputDir, connected)
% Subtracts median of each channel, then subtracts median of each time
% point.
%
% filename should include the extension
% outputDir is optional, by default will write to the directory of the input file
%
% should make chunk size as big as possible so that the medians of the
% channels differ little from chunk to chunk.
% 
% MGC 12/3/2018
% added an extra input which allows median of a subset of channels 
% (for when most of the probe is broken)
% connected: an extra input that specifies which channels to include in the
% time-point median step

chunkSize = 1000000;

fid = []; fidOut = [];

d = dir(filename);
nSampsTotal = d.bytes/nChansTotal/2;
nChunksTotal = ceil(nSampsTotal/chunkSize);
try
  
  [pathstr, name, ext] = fileparts(filename);
  fid = fopen(filename, 'r');
  if nargin < 3
    outputFilename  = [pathstr filesep name '_CAR' ext];
    mdTraceFilename = [pathstr filesep name '_medianTrace.mat'];
  else
    outputFilename  = [outputDir filesep name '_CAR' ext];
    mdTraceFilename = [outputDir filesep name '_medianTrace.mat'];
  end
  fidOut = fopen(outputFilename, 'w');
  
  % theseInds = 0;
  chunkInd = 1;
  medianTrace = zeros(1, nSampsTotal);
  while 1
    
    fprintf(1, 'chunk %d/%d\n', chunkInd, nChunksTotal);
    
    dat = fread(fid, [nChansTotal chunkSize], '*int16');
    
    if ~isempty(dat)
      
      %         theseInds = theseInds(end):theseInds(end)+chunkSize-1;
      
      dat = bsxfun(@minus, dat, median(dat,2)); % subtract median of each channel
      tm = median(dat(connected,:),1); % MALCOLM EDIT: only uses connected channels
      dat = bsxfun(@minus, dat, tm); % subtract median of each time point
      fwrite(fidOut, dat, 'int16');
      medianTrace((chunkInd-1)*chunkSize+1:(chunkInd-1)*chunkSize+numel(tm)) = tm;
      
    else
      break
    end
    
    chunkInd = chunkInd+1;
  end
  
  save(mdTraceFilename, 'medianTrace', '-v7.3');
  fclose(fid);
  fclose(fidOut);
  
catch me
  
  if ~isempty(fid)
    fclose(fid);
  end
  
  if ~isempty(fidOut)
    fclose(fidOut);
  end
  
  
  rethrow(me)
  
end
