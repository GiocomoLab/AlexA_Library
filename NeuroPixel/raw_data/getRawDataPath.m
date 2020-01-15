function [raw_path] = getRawDataPath(name,fileNameStruct)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for iF=1:numel(fileNameStruct)
    if endsWith(fileNameStruct(iF).name,name)
       raw_path = fullfile(fileNameStruct(iF).folder,fileNameStruct(iF).name);
       return
    end
end
raw_path = '';
end

