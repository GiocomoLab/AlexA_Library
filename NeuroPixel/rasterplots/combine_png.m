% script to combine rasters into one giant png for each session
% MGC 3/1/2019

%% params

% size of images

numrow = 8; % number of rows in final image

% session names
sn = dir('Z:\giocomo\attialex\images\rastersRegionParent_subset20');
session_name = {};
for iS=1:numel(sn)
    if sn(iS).isdir && ~ismember('.',sn(iS).name)
    session_name{end+1}=sn(iS).name;
    end
end
% where to save images
image_save_dir = 'F:\images\pretty_rasters_whole_session_parentSubset20';
mkdir(image_save_dir)
%% iterate over sessions
for k = 1:numel(session_name)
    
    %image_dir = fullfile('F:\images\',session_name{k},'\pretty_rasters\');
    image_dir = fullfile('Z:\giocomo\attialex\images\rastersRegionParent_subset20',session_name{k});

    % get png file names
    png_files = dir(fullfile(image_dir,'*.png'));
    png_files = {png_files.name};
    if numel(png_files)>300
        iidx = randsample(1:numel(png_files),300);
        png_files = png_files(iidx);
    end

    % find the files that are of the form (number).png
    %contains_number = regexp(png_files,'(\d*).png');
   % png_files = png_files(cell2mat(contains_number)==1);

%     % sort by cell number
%     cell_num = nan(numel(png_files),1);
%     for i = 1:numel(png_files)
%         this_split = strsplit(png_files{i},'.');
%         cell_num(i) = str2num(this_split{1});
%     end
%     [cell_num,sort_idx] = sort(cell_num);
%     png_files = png_files(sort_idx);
    dat = imread(fullfile(image_dir,png_files{1}),'png');
    [size_vert,size_horiz,~]=size(dat);


    % create empty matrix for holding final image
    numcol = ceil(numel(png_files)/numrow);
    final_image = nan(numrow*size_vert,numcol*size_horiz,3);

    % fill in final_image with data from png files
    for i = 1:numel(png_files)
        fprintf('session %d/%d: %s, file %d/%d\n',k,numel(session_name),session_name{k},i,numel(png_files));

        % read image
        dat = imread(fullfile(image_dir,png_files{i}),'png');

        % get row number and column number
        row = ceil(i/numcol);
        col = mod(i-1,numcol)+1;

        % enter data into final image
        final_image((row-1)*size_vert+1:row*size_vert,(col-1)*size_horiz+1:col*size_horiz,:) = dat;
    end
    
    % convert to uint8
    final_image = uint8(final_image);
    
    % write to file
    imwrite(final_image,fullfile(image_save_dir,sprintf('%s_all_rasters_combined.png',session_name{k})));
end