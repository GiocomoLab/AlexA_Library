function [data,shift]=correct_line_shift(data,template,varargin)
%%[data,shift]=correct_line_shift(data,template,varargin)
% estimates and corrects the shift between successive scan lines
%
%   data:          image stacks loaded in matrix format.
% optional parameters:
%   template:      mean of data if template is not specified.
%   'shift',value: specified shift value

p = inputParser;
if verLessThan('matlab','8.2') % not tested
    p.addParamValue('shift',[],@isnumeric)
    p.parse(varargin{:})
else
    addRequired(p,'data',@isnumeric)
    addOptional(p,'template',[],@isnumeric)
    addParameter(p,'shift',[],@isnumeric)
    if exist('template','var')
        parse(p,data,template,varargin{:})
    else
        parse(p,data,varargin{:})
    end
end
data = p.Results.data;
template = p.Results.template;
shift = p.Results.shift;


if isempty(template) && isempty(shift)
    template=mean(data,3);
end

twoddata=0;
if size(data,3)==1
    tdata(:,:,1)=data;
    tdata(:,:,2)=data;
    twoddata=1;
    data=tdata;
end

if isempty(shift)
    shift=estimate_line_shift(template);
    %disp(['Estimated shift is ' num2str(shift)]);
end

if shift~=0
    for ind=1:size(data,1)
        if rem(ind,2)==1
            data(ind,:,:)=circshift(squeeze(data(ind,:,:)),[shift 0]);
        end
    end
end

if twoddata
    data=data(:,:,1);
end

