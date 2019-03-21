% directoryPath = 'F:\H1_flipped\2019-03-12_10-21-28'
function neuralynx2kilosort(directoryPath,target_dir)
    % takes a directory path (with no trailing backslash) and outputs a
    % .bin file that is a data matrix with samples of size 64channels x numOfSamples
    if nargin ==1
        target_dir = directoryPath;
    end
    
    fprintf(strcat('\nStart Processing: ', datestr(now,'mmmm dd, yyyy HH:MM:SS AM'),'\n'));
    numOfChannels = 64;
    dataMatrix = zeros(numOfChannels,100,'int16');
    for csc = 1:numOfChannels
        cscPath = fullfile(directoryPath, ['CSC_HP_' num2str(csc) '_0001.ncs']);
        %cscPath = fullfile(directoryPath, ['CSC' num2str(csc) '.ncs']);
        [Samples,header]=Nlx2MatCSC(cscPath, [0 0 0 0 1], 1, 1, [] );
        tmp=split(header{17}); %assuming conversion factor is in here;
        conv_factor = str2double(tmp{2}); % in volts
        conv_factor = conv_factor*10e6; % in micro volts
        %conv_factor = 1;
        Samples = reshape(Samples,1,[])*-1*conv_factor; 
        
        
         Samples= int16(Samples); %load neuralynx file, linearize samples, convert to int16
        if csc == 1 %pre-allocate dataMatrix size on first pass to increase speed
           numOfSamples = size(Samples,2);
           dataMatrix = zeros(numOfChannels,numOfSamples,'int16');
        end
        dataMatrix(csc,:) = Samples; %allocate csc data to appropriate row in data matrix
        fprintf(strcat('\nProcessed: ',num2str(csc), ' out of 64 CSC files.'));
    end
    [~,Name,~] = fileparts(directoryPath);
    fid = fopen(fullfile(target_dir,[Name '_dataMatrix.bin']), 'w'); 
    fwrite(fid, dataMatrix, 'int16');
    fclose(fid);
    fprintf(strcat('\nDone Processing: ', datestr(now,'mmmm dd, yyyy HH:MM:SS AM'),'\n'));
end


