function sbx2tif(fname,varargin)

% sbx2tif
% Generates tif file from sbx files
% Argument is the number of frames to convert
% If no argument is passed the whole file is converted
% edited by Mai-Anh Vu, 1/3/2019, to automatically detect whether multiple
% channels were recorded and to save out tif files accordingly

z = sbxread(fname,1,1);
global info;

if(nargin>1)
    N = min(varargin{1},info.max_idx);
else
    N = info.max_idx;
end
nchan = info.nchan;

for chan = 1:nchan
    k = 0;
    done = 0;
    while(~done && k<=N)
        try
            q = sbxread(fname,k,1);
            q = squeeze(q(chan,:,:));
            if(k==1)
                imwrite(q,[fname '_chan' num2str(chan) '.tif'],'tif');
            else
                imwrite(q,[fname '_chan' num2str(chan) '.tif'],'tif','writemode','append');
            end
        catch
            done = 1;
        end
        k = k+1;
    end
end