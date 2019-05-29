function final_image = combine_images(image_dir,n_rows)

% script to combine rasters into one giant png for each session
% MGC 3/1/2019

%% params

% size of images

numrow = n_rows; % number of rows in final image



   

    % get png file names
    png_files = dir(sprintf('%s\\*.png',image_dir));
    png_files = {png_files.name};

    dat = imread(fullfile(image_dir,png_files{1}),'png');
    size_vert = size(dat,1);
    size_horiz = size(dat,2);

    % sort by cell number
    cell_num = nan(numel(png_files),1);
    cell_depth = nan(numel(png_files),1);
    for ic = 1:numel(png_files)
        this_split = strsplit(png_files{ic},{'.','_'});
        
        cell_num(ic) = str2double(this_split{1});
        cell_depth(ic) = str2double(this_split{2});
        
    end
    [cell_num,sort_idx] = sort(cell_depth);
    png_files = png_files(sort_idx);

    % create empty matrix for holding final image
    numcol = ceil(numel(png_files)/numrow);
    final_image = nan(numrow*size_vert,numcol*size_horiz,3);

    % fill in final_image with data from png files
    for i = 1:numel(png_files)
        %fprintf('session %d/%d: %s, file %d/%d\n',k,numel(session_name),session_name{k},i,numel(png_files));

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

end

