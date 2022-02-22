 function view_stack(data, varargin)
%%view_stack(data, varargin)
% viewer for stacks and movies
%
%   data - our stack
%
% optional settings (varargin)
%   map - optional change of colormap, default gray
%   adata - aux data to plot on the right of the image (e.g. running speed)
%   FrameRate - 0 (default) or FrameRate in Hz (e.g. 15 for 4 pizo layers on the
%      12khz setup - display ms instead of frame number
%
% key control:
% shift+arrow up/down - adjust brightness
% arrow up/down - adjust movie speed
% space tab: pause movie
% a: show activity of selected area
% f: mark line in figure
% c: clear all ROIs
% r: compare first and last averaged frames
% s: scale to square
% m: make movie
% t: show average image of stack
% v: sliding average over n frames
% w: add mark in stack (e.g. to label active cells)
% W: delete mark in stack
% x: toggle export frame parameter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p=inputParser;

addRequired(p,'data',@isnumeric)

    addParameter(p,'marks',[],@isnumeric)
    addParameter(p,'adata',[],@isnumeric)
    addParameter(p,'FrameRate',0,@isnumeric)
    addParameter(p,'FrameZero',0,@isnumeric)
    addParameter(p,'map','gray')
    addParameter(p,'targetsXY',[],@isnumeric)
    addParameter(p,'clim',[])



parse(p,data,varargin{:})

data=p.Results.data;
params.marks=p.Results.marks;
map=p.Results.map;
params.adata=p.Results.adata;
params.FrameRate=p.Results.FrameRate;
params.FrameZero=p.Results.FrameZero;
params.clim = p.Results.clim;

if isempty(params.adata)
    params.draw_adata = 0;
else
    params.draw_adata = 1;
    params.adata(params.adata<0) = 0;
    params.adata = ntzo(params.adata);
end

if ~isempty(p.Results.targetsXY)
    drawtargets = true;
    targetsXY = p.Results.targetsXY;
else
    drawtargets = false;
end
 
hf=figure;
clf;
params.button_down=0;
params.movie_status=0;
params.frame_pause=0.01; % pause in seconds between frames in playback
params.playback_spacing=1; % only every n-th frame is shown during playback default=1 (every frame)
params.as_ind=0; % index over area-selected for plotting activity
params.nbr_avg_frames=1;
params.tx=0; % tilt of image in x
params.ty=0; % tilt of image in y
params.marksh=[];
params.export_frame=0;

set(hf,'UserData',params);
dimensions = size(data);

offset = 200;

set(hf,'Position',[offset offset dimensions(2) 1.03*dimensions(1)]);

frame_ind=1;

template=data(:,:,round(size(data,3)/2));

low_contrast_lim=prctile(reshape(template,numel(template),1),10);
high_contrast_lim=prctile(reshape(template,numel(template),1),99.99);

if isempty(params.clim)
    if low_contrast_lim< high_contrast_lim
        params.clim=[low_contrast_lim high_contrast_lim];
    else
        disp('set params.clim to default values')
        params.clim=[0 255];
    end
end
if params.draw_adata
    lim1 = 0.975;
    params.h_aux_ax = axes('position',[0.975 0.03 0.025 0.97],'color','k');
    params.h_aux_data = plot([0 0],[0 params.adata(1)],'linewidth',10,'color','g');
    set(params.h_aux_ax,'Visible','off','YLim',[0 1])
else
    lim1 = 1;
end
params.h_f_ax=axes('position',[0 0.03 lim1 0.97],'color','k','clim',params.clim);
params.h_im_data=imagesc(data(:,:,frame_ind));
set(params.h_f_ax,'Visible','off')

hold on
for ind=1:length(params.marks)
    params.marksh(end+1)=plot(-100,-100,'r.','linewidth',2,'markersize',10);
end

% set(params.h_im_data,'CDataMapping','direct');

yl=ylim;xl=xlim;
params.h_txt=text(xl(2)/20,yl(2)/20,'1','fontsize',12,'color','w','fontweight','bold');

colormap(map);
if drawtargets
    for t = 1:size(targetsXY,1)
        if targetsXY(t,3)
            tColor = 'r';
            markerform = 'o';
        else
            tColor = 'c';
            markerform = '-s';
        end
        plot(targetsXY(t,1),targetsXY(t,2),markerform,'MarkerSize',6,'color',tColor)
    end
end

h_sl_ax=axes('position',[0 0 1 0.03],'color','k');
sl_x_pos=0;
h_sl_bar=plot([1 1]*sl_x_pos,[0 1],'r','linewidth',20);
set(h_sl_ax,'Visible','off','XLim',[0 1],'YLim',[0 1])

set(hf,'windowbuttondownfcn',@view_stack_buttdofcn);
set(hf,'windowbuttonupfcn',@view_stack_buttupfcn);
set(hf,'windowbuttonmotionfcn',{@view_stack_winmotfcn,h_sl_ax,h_sl_bar,data});
set(hf,'KeyPressFcn' ,{@view_stack_keypress,h_sl_bar,data});
set(hf,'CloseRequestFcn',@close_fcn);
set(hf,'menubar','none','color','k');
set(params.h_f_ax,'clim',params.clim);
set(params.h_f_ax,'TickLength',[0 0]);
set(hf,'UserData',params);

function view_stack_winmotfcn(hf,e,h_sl_ax,h_sl_bar,data)
params=get(hf,'UserData');
if params.button_down
    cp=get(h_sl_ax,'currentpoint');
    cp=cp(1);
    cp=min(cp,1);
    cp=max(cp,0);
    draw_frame(cp,h_sl_bar,params,data);
end

function view_stack_buttdofcn(hf,e)
params=get(hf,'UserData');
params.button_down=1;
set(hf,'UserData',params);

function view_stack_buttupfcn(hf,e)
params=get(hf,'UserData');
params.button_down=0;
set(hf,'UserData',params);

function view_stack_keypress(hf,event,h_sl_bar,data)
params=get(hf,'UserData');
nFrames = size(data,3);
spacing = 1/nFrames;
switch event.Character
    case 's'  % scale
        temp = get(hf);
        figLength = temp.Position(3);
        set(hf,'Position',[temp.Position(1) temp.Position(2) figLength 1.03*figLength*size(data,1)/size(data,2)]);
        
        
    case 'h' % help
         help view_stack

    case 'm' % make a movie
        [avi_fname,avi_path]=uiputfile('tiff_stack.avi','save stack as');
        mov=VideoWriter([avi_path avi_fname],'Motion JPEG AVI');

        mov.FrameRate = 23;
        speed_factor=input('Select speed to save movie at: ');
        frame_bounds=input('Select [start_frame stop_fram], 0 for all frames: ');
        open(mov);
        if frame_bounds==0
            frame_bounds=[0 nFrames];
        end
        cp2 = ntzo(frame_bounds(1)/nFrames:spacing*speed_factor:(frame_bounds(2)+1)/nFrames-spacing*speed_factor);
        ii = 1;
        for cp=frame_bounds(1)/nFrames:spacing*speed_factor:(frame_bounds(2)+1)/nFrames-spacing*speed_factor
            
            curr_frame_ind=round(cp*(nFrames-1))+1;
            set(h_sl_bar,'Xdata',[1 1]*cp2(ii));
            
            set(params.h_im_data,'CData',mean(data(:,:,curr_frame_ind:min(size(data,3),curr_frame_ind-1+speed_factor+params.nbr_avg_frames)),3));
            
            if params.FrameRate==0
                set(params.h_txt,'string',num2str(curr_frame_ind));
            else
                set(params.h_txt,'string',[num2str(floor((curr_frame_ind-params.FrameZero)/params.FrameRate)) ' s']);
            end
            
            if params.draw_adata
                set(params.h_aux_data,'YData',[0 params.adata(curr_frame_ind)]);
            end
            frame = getframe(gcf);
            writeVideo(mov,frame);
            ii = ii + 1;
        end
        close(mov);
    case 'a' % show activity of selected area
        params.as_ind=params.as_ind+1;
        color_ind='mgycrbw';
        [asx,asy]=ginput(2);%,'Color',[1 1 1]);
        hold on
%         plot([as(2) as(2) as(1) as(1) as(2)],[as(4) as(3) as(3) as(4) as(4)],color_ind(params.as_ind));
        plot([asx(2) asx(2) asx(1) asx(1) asx(2)],[asy(2) asy(1) asy(1) asy(2) asy(2)],color_ind(mod(params.as_ind-1,length(color_ind))+1));
%         as_x=round(sort([as(3) as(4)]));
%         as_y=round(sort([as(1) as(2)]));
        as_x=round(sort([asy(1) asy(2)]));
        as_y=round(sort([asx(1) asx(2)]));
        fig_pos=get(gcf,'position');
        if isfield(params,'sup_fig_h') && ishandle(params.sup_fig_h)
            set(params.sup_fig_h,'position',[fig_pos(1) fig_pos(2)-100-33 fig_pos(3) 100]);
        else
            params.sup_fig_h=figure('position',[fig_pos(1) fig_pos(2)-100-33 fig_pos(3) 100],'menubar','none','color','k');
            axes('position',[0 0 1 1],'color','k');
            hold on
        end
        prev_fig_handle=gcf;
        figure(params.sup_fig_h);
        try
            raw_act_trace=squeeze(mean(mean(data(as_x(1):as_x(2),as_y(1):as_y(2),:),2),1));
            plot(raw_act_trace,color_ind(mod(params.as_ind-1,length(color_ind))+1));
            assignin('base','raw_act_trace',raw_act_trace);
            evalin('base',['raw_traces(:,' num2str(params.as_ind) ')=raw_act_trace;']);
        catch
            disp('try again - selection not valid');
        end
        axis tight;
        figure(prev_fig_handle);
        set(hf,'UserData',params);
    case 'f' % mark line in figure
        as=ginput(2);
        hold on
        plot([as(1) as(2)],[as(3) as(4)],'r','linewidth',2)
        
    case 't' % show template
        fig_pos=get(gcf,'position');
        prev_fig_handle=gcf;
        if isfield(params,'sup_fig2_h') && ishandle(params.sup_fig2_h)
            set(params.sup_fig2_h,'position',[fig_pos(1)+fig_pos(3) fig_pos(2) fig_pos(3) fig_pos(4)]);
            figure(params.sup_fig2_h);
        else
            params.sup_fig2_h=figure('position',[fig_pos(1)+fig_pos(3) fig_pos(2) fig_pos(3) fig_pos(4)],'menubar','none');
        end
        axes('position',[0 0 1 1]);
        
        imagesc(mean(data,3));
        
        colormap gray;
        %         set(gca,'clim',params.clim)
        axis off
        figure(prev_fig_handle);
        set(hf,'UserData',params);
    case 'c' % clear
        hold off
        cp=get(h_sl_bar,'Xdata');
        cp=cp(1);
        axes(params.h_f_ax);
        params.h_im_data=imagesc(data(:,:,round(cp*(size(data,3)-1))+1));
        set(params.h_f_ax,'clim',params.clim,'xtick',[]);
        yl=ylim;xl=xlim;
        params.h_txt=text(xl(2)/20,yl(2)/20,num2str(round(cp*(size(data,3)-1))+1),'fontsize',12,'color','w','fontweight','bold');
        if isfield(params,'sup_fig_h') && ishandle(params.sup_fig_h)
            close(params.sup_fig_h);
        end
        if isfield(params,'sup_fig2_h') && ishandle(params.sup_fig2_h)
            close(params.sup_fig2_h);
        end
        params.as_ind=0;
        set(hf,'UserData',params);
        draw_frame(cp,h_sl_bar,params,data);
        evalin('base',['clear raw_traces;']);
    case 'v' % set number of frames to average
        params.nbr_avg_frames=input('How many frames do you want to average: ');
        set(hf,'UserData',params);
        figure(hf);
    case 'x' % toggle export frame parameter
        params.export_frame=~params.export_frame;
        if params.export_frame==1
            disp('Export frames is turned ON');
        else
            disp('Export frames is turned OFF');
        end
        set(hf,'UserData',params);
    case 'r' % compare first and last frames
        params.nbr_insp_frames=input('How many frames to use for inspection: ');
        set(hf,'UserData',params);
        fig_pos=get(gcf,'position');
        prev_fig_handle=gcf;
        if isfield(params,'sup_fig2_h') && ishandle(params.sup_fig2_h)
            set(params.sup_fig2_h,'position',[fig_pos(1)+fig_pos(3) fig_pos(2) fig_pos(3) fig_pos(4)]);
            figure(params.sup_fig2_h);
        else
            params.sup_fig2_h=figure('position',[fig_pos(1)+fig_pos(3) fig_pos(2) fig_pos(3) fig_pos(4)],'menubar','none');
        end
        axes('position',[0 0 1 1]);
        %         [shift_x,shift_y] = register_frames(data(:,:,[1:params.nbr_insp_frames end-params.nbr_insp_frames:end]),mean(data(:,:,10:20),3),0.2);
        [shift_x,shift_y] = register_frames(data(:,:,[1:params.nbr_insp_frames 100:100+params.nbr_insp_frames]),mean(data(:,:,10:20),3),0.2);
        
        tempdata=shift_data(data(:,:,[1:params.nbr_insp_frames end-params.nbr_insp_frames:end]),shift_x,shift_y);
        %         tempdata = data(:,:,[1:params.nbr_insp_frames end-params.nbr_insp_frames:end]);
        imshowpair(mean(tempdata(:,:,1:params.nbr_insp_frames),3), mean(tempdata(:,:,params.nbr_insp_frames+1:end),3),'montage')
        set(gca,'clim',params.clim)
        axis off
        figure(prev_fig_handle);
        set(hf,'UserData',params);
    case 'w' % add mark in stack (e.g. to label active cells)
        mark=round(ginput(1));
        cp=get(h_sl_bar,'Xdata');
        cp=cp(1);
        cf=round(cp*(size(data,3)-1))+1;
        params.marks(end+1)=size(data,1)*size(data,2)*(cf-1)+(mark(1)-1)*size(data,1)+mark(2);
        hold on
        params.marksh(end+1)=plot(mark(1),mark(2),'ro','linewidth',2,'markersize',8);
        set(hf,'UserData',params);
        assignin('base','marks',params.marks)
    case 'W' % delete mark in stack
        mark=round(ginput(1));
        cp=get(h_sl_bar,'Xdata');
        cp=cp(1);
        cf=round(cp*(size(data,3)-1))+1;
        for ind=1:length(params.marks)
            [tmpx,tmpy,tmpz]=ind2sub(size(data),params.marks(ind));
            if abs(tmpx-mark(2))<3 && abs(tmpy-mark(1))<3 && abs(tmpz-cf)<3
                delete(params.marksh(ind));
                params.marksh(ind)=[];
                params.marks(ind)=[];
                break
            end
        end
        set(hf,'UserData',params);
        assignin('base','marks',params.marks)
        
    case 'p'
        assignin('base','movie_params',params);
end

switch event.Key
    case 'rightarrow'
        if strcmp(event.Modifier,'shift')
            if params.clim(1)==0
                params.clim(1)=1;
            end
            params.clim(1)=params.clim(1)*1.1;
            set(params.h_f_ax,'clim',params.clim)
            set(hf,'UserData',params);
        elseif strcmp(event.Modifier,'control')
            params.ty=params.ty+1;
            set(hf,'UserData',params);
            cp=get(h_sl_bar,'Xdata');
            cp=cp(1);
            draw_frame(cp,h_sl_bar,params,data);
        else
            cp=get(h_sl_bar,'Xdata');
            cp=cp(1);
            cp = cp + spacing;
            if cp > 1, cp = 1; end;
            draw_frame(cp,h_sl_bar,params,data);
        end
    case 'leftarrow'
        if strcmp(event.Modifier,'shift')
            params.clim(1)=params.clim(1)*0.9;
            set(params.h_f_ax,'clim',params.clim)
            set(hf,'UserData',params);
        elseif strcmp(event.Modifier,'control')
            params.ty=params.ty-1;
            set(hf,'UserData',params);
            cp=get(h_sl_bar,'Xdata');
            cp=cp(1);
            draw_frame(cp,h_sl_bar,params,data);
        else
            cp=get(h_sl_bar,'Xdata');
            cp=cp(1);
            cp = cp - spacing;
            if cp < 0, cp = 0; end;
            draw_frame(cp,h_sl_bar,params,data);
        end
    case 'uparrow'
        if strcmp(event.Modifier,'shift')
            params.clim(2)=params.clim(2)*1.1;
            set(params.h_f_ax,'clim',params.clim)
        elseif strcmp(event.Modifier,'control')
            params.tx=params.tx+1;
            set(hf,'UserData',params);
            cp=get(h_sl_bar,'Xdata');
            cp=cp(1);
            draw_frame(cp,h_sl_bar,params,data);
        else
            params=get(hf,'UserData');
            if params.frame_pause<=0.01
                params.playback_spacing=params.playback_spacing*2;
            else
                params.frame_pause=params.frame_pause*0.5;
            end
        end
        set(hf,'UserData',params);
        
    case 'downarrow'
        if strcmp(event.Modifier,'shift')
            params.clim(2)=params.clim(2)*0.9;
            set(params.h_f_ax,'clim',params.clim)
        elseif strcmp(event.Modifier,'control')
            params.tx=params.tx-1;
            set(hf,'UserData',params);
            cp=get(h_sl_bar,'Xdata');
            cp=cp(1);
            draw_frame(cp,h_sl_bar,params,data);
        else
            params=get(hf,'UserData');
            if params.playback_spacing>1
                params.playback_spacing=params.playback_spacing/2;
            else
                params.frame_pause=params.frame_pause*2;
            end
        end
        set(hf,'UserData',params);
    case 'space'   % play as movie
        params.movie_status=~params.movie_status;
        set(hf,'UserData',params);
        cp=get(h_sl_bar,'Xdata');
        cp=cp(1);
        while params.movie_status
            params=get(hf,'UserData');
            cp = cp + spacing*params.playback_spacing;
            if cp > 1
                cp = 0;
            end;
            draw_frame(cp,h_sl_bar,params,data);
            drawnow
            pause(params.frame_pause);
        end
end

function draw_frame(cp,h_sl_bar,params,data)
curr_frame_ind=round(cp*(size(data,3)-1))+1;
set(h_sl_bar,'Xdata',[1 1]*cp);

if params.nbr_avg_frames==1
    if params.tx~=0 || params.ty~=0
        [pix_inds]=get_rot_ind(size(data),curr_frame_ind,params.tx,params.ty);
        frame=reshape(data(pix_inds),[size(data,1) size(data,2)]);
    else
        frame=data(:,:,curr_frame_ind);
    end
    set(params.h_im_data,'CData',frame);
    for ind=1:length(params.marks)
        if ~isempty(intersect(pix_inds,[-2:2]*(size(data,1)*size(data,2))+params.marks(ind)))
            [tmpx,tmpy,~]=ind2sub(size(data),params.marks(ind));
            set(params.marksh(ind),'XData',tmpy,'YData',tmpx);
        else
            set(params.marksh(ind),'XData',-100,'YData',-100);
        end
    end
else
    if params.tx~=0 || params.ty~=0
        frames=zeros(size(data,1),size(data,2),params.nbr_avg_frames);
        for ind=1:params.nbr_avg_frames
            [pix_inds]=get_rot_ind(size(data),curr_frame_ind-1+ind,params.tx,params.ty);
            frames(:,:,ind)=reshape(data(pix_inds),[size(data,1) size(data,2)]);
        end
        frame=max(frames,[],3);
    else
        frame=mean(data(:,:,curr_frame_ind:min(curr_frame_ind+params.nbr_avg_frames-1,size(data,3))),3);
    end
    set(params.h_im_data,'CData',frame);
    
end
if params.draw_adata
    set(params.h_aux_data,'YData',[0 params.adata(curr_frame_ind)]);
end
if params.export_frame
    assignin('base','curr_frame',frame);
end

if params.tx==0 && params.ty==0
    if params.FrameRate==0
        set(params.h_txt,'string',num2str(curr_frame_ind));
    else
        set(params.h_txt,'string',[num2str(floor((curr_frame_ind-params.FrameZero)/params.FrameRate)) ' s']);
    end
else
    set(params.h_txt,'string',[num2str(curr_frame_ind) ' - Tilt X: ' num2str(params.tx) ' Tilt Y: ' num2str(params.ty)]);
end

function close_fcn(hf,e)
params=get(hf,'UserData');
if isfield(params,'sup_fig_h') && ishandle(params.sup_fig_h)
    close(params.sup_fig_h);
end
if isfield(params,'sup_fig2_h') && ishandle(params.sup_fig2_h)
    close(params.sup_fig2_h);
end
delete(hf)