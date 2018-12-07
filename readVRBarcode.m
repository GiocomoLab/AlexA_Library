%% Import data from text file.
% Script for importing data from the following text file:
%
%    X:\Projects\MismatchExperiment\E2\VR\0827_mismatch_np_1_position.txt
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2018/08/31 11:51:48

%% Initialize variables.
filename = 'F:\G4\1204_gaincontrast_2\1204_gaincontrast_2_position.txt';
delimiter = '\t';

%% Format for each line of text:
%   column5: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
timestamp = dataArray{:,4};
barcode1 = dataArray{:, 5};


%% Clear temporary variables
clearvars filename delimiter formatSpec fileID ans;