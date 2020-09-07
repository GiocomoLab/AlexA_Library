function ca(listOfFigsToClose)
% close all figures - but not file references, GUIs etc.
% listOfFigsToClose can be the range of figures (i.e. can contain values
% that are not figures) to close
% useage: ca(1:4) - closes figures 1 through 4. does not return an error if
% e.g. 3 does not exist.

if nargin<1
    listOfFigsToClose=[];
end

figHandles = get(0,'Children');
figHandles = findobj('Type','figure');


if isempty(listOfFigsToClose)
    for ind=1:length(figHandles)
        if figHandles(ind)~=1001
            try
                % view stack closes subfigures automatically on close of the
                % main figure
                close(figHandles(ind));
            end
        end
    end
else
    for ind=1:length(figHandles)
        if sum(figHandles(ind)==listOfFigsToClose)>0
            try
                % view stack closes subfigures automatically on close of the
                % main figure
                close(figHandles(ind));
            end
        end
    end
end