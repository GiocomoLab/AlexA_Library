function app=POINav_createComponents(app,nRows,nCols)
if nargin ==1
    nRows=256;
    nCols=256;
end
% Create Figure
app.hf.Position = [100 100 1000 680];

% create userdata variables
UserData=app.hf.UserData;
UserData.OVERLAY=0;
UserData.selectedRow=1;
UserData.N_Images= 15;
UserData.templates=zeros(nRows,nCols,UserData.N_Images);
UserData.commentColumn=7;

app.hf.UserData=UserData;

% Create Live_Ax
app.Live_Ax = axes(app.hf);
app.Live_Ax.Units='pixels';
title(app.Live_Ax, 'Live');
%app.Live_Ax.Box = 'on';
%app.Live_Ax.XColor = [0.149 0.149 0.149];
%app.Live_Ax.XColorMode = 'manual';
app.Live_Ax.Position = [10 385 400 269];
imageData=randn(nRows,nCols);
app.Live_Ax_Image=imagesc(app.Live_Ax,imageData);
set(app.Live_Ax,'xtick',[],'ytick',[]);


% Create Template_Ax
app.Template_Ax = axes(app.hf);
app.Template_Ax.Units='pixels';

title(app.Template_Ax, 'Template');

app.Template_Ax.Position = [409 385 400 269];
imageData=zeros(nRows,nCols);
app.Template_Ax_Image=imagesc(app.Template_Ax,imageData);
set(app.Template_Ax,'xtick',[],'ytick',[]);

% Create coord system_ax
app.Map_Ax = axes(app.hf);
app.Map_Ax.Units='pixels';

title(app.Map_Ax, 'Map');

app.Map_Ax.Position = [700 80 250 250];
%imageData=zeros(nRows,nCols);
%app.Template_Ax_Image=imagesc(app.Template_Ax,imageData);
%set(app.Template_Ax,'xtick',[],'ytick',[]);

% Create UITable
app.UITable = uitable(app.hf);
app.UITable.ColumnName = {'X'; 'Y'; 'Z';'dx';'dy';'dz'; 'Comment'};
%app.UITable.RowName = {};
app.UITable.Position = [58 144 613 185];
d = {0,0,0,0,0,0,''};

f=repmat(d,10,1);
app.UITable.Data=f;

app.UITable.CellEditCallback={@TableEditCallback,app};
app.UITable.CellSelectionCallback = {@TableSelectCallback,app};
app.UITable.ColumnEditable=true(1,7);

app.FTA = uicontrol('style','edit','string',int2str(app.hf.UserData.FTA),...
    'callback',{@FTASet,app});
app.FTA.Position = [110 348 30 22];
app.FTALabel=uicontrol('style','text','string','FTA');
app.FTALabel.Position=[80 348 30 22];
% Create SetPOIButton
app.SetPOIButton = uicontrol('style','push',...
    'string','Set POI',...
    'callback',{@SetPOIButtonPushed, app});
app.SetPOIButton.Position = [140 348 100 22];

% Create SetORIButton

app.SetPOIButton = uicontrol('style','push',...
    'string','Set ORI',...
    'callback',{@SetORIButtonPushed, app});
app.SetPOIButton.Position = [281 348 100 22];

% Create SetLoadButton

app.SetPOIButton = uicontrol('style','push',...
    'string','Load',...
    'callback',{@SetLoadButtonPushed, app});
app.SetPOIButton.Position = [414 348 100 22];

% Create SetLoadButton

app.SetPOIButton = uicontrol('style','push',...
    'string','Save',...
    'callback',{@SetSaveButtonPushed, app});
app.SetPOIButton.Position = [559 348 100 22];

% Create toggleButton
app.OverlayButton = uicontrol('style','checkbox',...
    'units','pixels',...
    'position',[659 348 100 22],...
    'string','Merge',...
    'val',app.hf.UserData.OVERLAY,...
    'callback',{@OverlayButtonPushed, app});




end

function FTASet(varargin)
app=varargin{3};
fta=varargin{1}.String;
app.hf.UserData.FTA=str2num(fta);

end

function TableSelectCallback(varargin)
ed = varargin{2};
app = varargin{3};
if ~isempty(ed.Indices)
    app.hf.UserData.selectedRow=ed.Indices(1);
end
end

function plotCoordinates(app)
%ud=app.hf.UserData;
table=app.UITable.Data;
axes(app.Map_Ax);
hold on;
for ii=1:10
    x=table{ii,1};
    y=table{ii,2};
    
    if x~=0
        
        plot(x,y,'.')
        text(x,y,int2str(ii),'VerticalAlignment','bottom','HorizontalAlignment','right')
    end
end
end

% Button pushed function: SetPOIButton
function SetPOIButtonPushed(varargin)
app=varargin{3};
saveTMP(app)

ud=app.hf.UserData;
CData=app.Live_Ax_Image.CData;
x=ud.currX;
y=ud.currY;
z=ud.currZ;
selectedRow=app.hf.UserData.selectedRow;
%check if empty
%assignin('base','table',app.table.Data
comment=app.UITable.Data{selectedRow,ud.commentColumn};
XX=app.UITable.Data{selectedRow,1};
if selectedRow~=1
oriX=app.UITable.Data{1,1};
oriY=app.UITable.Data{1,2};
oriZ=app.UITable.Data{1,3};
else
    oriX=x;
    oriY=y;
    oriZ=z;
end
dx=x-oriX;
dy=y-oriY;
dz=z-oriZ;
if isempty(comment) && XX==0
    update=true;
else
    answer = questdlg(sprintf('Overwrite POI in Row %d?',selectedRow), ...
        'ROI overwrite', ...
        'Yes','No','No');
    % Handle response
    switch answer
        case 'Yes'
            update=true;

            
        case 'No'
            update=false;
    end
    
end
if update
    ud.templates(:,:,selectedRow)=CData;
    app.Template_Ax_Image.CData=CData;
    app.UITable.Data{selectedRow,1}=x;
    app.UITable.Data{selectedRow,2}=y;
    app.UITable.Data{selectedRow,3}=z;

    app.UITable.Data{selectedRow,4}=dx;
    app.UITable.Data{selectedRow,5}=dy;
    app.UITable.Data{selectedRow,6}=dz;
    app.hf.UserData.selectedRow=selectedRow+1;
else
                warning('Roi not updated')

end
plotCoordinates(app);
end
% Button pushed function: SetORIButton
function SetORIButtonPushed(varargin)
app=varargin{3};
saveTMP(app)

end

% Button pushed function: LoadButton
function LoadButtonPushed(varargin)
app=varargin{3};
saveTMP(app)

end

% Button pushed function: SaveButton
function SetSaveButtonPushed(varargin)
app=varargin{3};
[file,path] = uiputfile('*.mat','POIFile');
tableData=app.UITable.Data;
templates=app.hf.UserData.templates;
saveTMP(app)
save([path file],'tableData','templates');
end

function TableEditCallback(varargin)
app=varargin{3};
disp('cell edited')
saveTMP(app)
end

function OverlayButtonPushed(varargin)
app=varargin{3};
ud=app.hf.UserData;
if ud.OVERLAY == 0
    ud.OVERLAY = 1;
else
    ud.OVERLAY = 0;
end
ud.OVERLAY
app.hf.UserData=ud;
%get(app.OverlayButton);
%set(handle,'val',app.OVERLAY);
end

function saveTMP(app)
fn=datestr(now,'mm-dd-yyyy_HH-MM-SS');
tableData=app.UITable.Data;
templates=app.hf.UserData.templates;
save(['C:\\temp\POINavigator\' fn '.mat'],'tableData','templates')
end
