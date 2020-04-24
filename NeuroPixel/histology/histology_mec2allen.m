histology = process_histology();
annotation_volume_location = 'F:\code\allenCCF\Allen\annotation_volume_10um_by_index.npy';
structure_tree_location = 'F:\code\allenCCF\Allen\structure_tree_safe_2017.csv';
template_volume_location = 'F:\code\allenCCF\Allen\template_volume_10um.npy';
afs = dir('Z:\giocomo\export\data\Projects\AlexA_NP\Histology\AA*');
animals = {afs.name};





%% GET AND PLOT PROBE VECTOR IN ATLAS SPACE

% load the reference brain annotations
if ~exist('av','var') || ~exist('st','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
end

%% plot in MEC coordinates
figure
hold on
for iP= 1:numel(histology.Mouse)
    
    pos3D = histology.origin(iP,:) + linspace(0,histology.FinalDepth(iP),10)'*histology.unit_vector(iP,:);
    probe_track = [pos3D(1,:);pos3D(end,:)];
    scatter3(pos3D(:,1),pos3D(:,2),pos3D(:,3))
    plot3(probe_track(:,1),probe_track(:,2),probe_track(:,3))
end
xlabel('ML')
ylabel('Layer')
zlabel('Depth')
axis equal

%% from MEC corrdinates to Allen Coordinates in MM:

fwireframe = plotBrainGrid([], [], [], true); hold on; 
fwireframe.InvertHardcopy = 'off';
brainfig = gcf;
probefig=figure();
set(gca,'ZDir','reverse')
hold on
controlfig = figure();
hold on
atlasRes = 0.010; % mm
bregma=allenCCFbregma()
MEC_Offset = [3.0,2.0,-5.25]; %ML,DV,AP %estimate from old paxinos atlas
MEC_angle = 15;
MEC_Offset = [3.3,2.0,-5.25]; % %ML,DV,AP from 2008 allen
AllPointsMEC = {};
for iP= 1:numel(histology.X1)
    mec_entry = [histology.X1(iP),histology.Y1(iP),histology.Z1(iP)]/1000;
    probe_term = [histology.X2(iP),histology.Y2(iP),histology.Z2(iP)]/1000;
    if any(isnan(probe_term))
        continue
    end
    entry_allen = mec_entry([1 3 2]) + MEC_Offset;
    t1=probe_term(1);
    t2=probe_term(3)*cos(deg2rad(20)); %play with this angle
    t3=probe_term(2)*sin(deg2rad(20))+probe_term(2);
    
    term_allen = [t1,t2,t3] + MEC_Offset; %ML,DV,AP
    pos3D=[entry_allen;term_allen]; %ML,DV,AP
    vecAllen = -term_allen + entry_allen;
    vecAllen = vecAllen/norm(vecAllen);
    if strcmp(histology.Side{iP},'L')
        %continue
        term_allen(1)=-term_allen(1); %flip to minus for left hemisphere
        pos3D(:,1)=-pos3D(:,1);
        vecAllen(1)=-vecAllen(1);
    end
    wholeTrack = [term_allen;term_allen+histology.FinalDepth(iP)*vecAllen/1000]; %ML,DV,AP
    figure(probefig)
    scatter3(pos3D(:,1),pos3D(:,3),pos3D(:,2),4,[1 0 0; 0 0 1]) %ML,AP,DV
    plot3(wholeTrack(:,1),wholeTrack(:,3),wholeTrack(:,2),'LineWidth',2,'Color','b') % plot x: ml, y: Ap, z:dv
    %plot3(probe_track(:,1),probe_track(:,3),probe_track(:,2))
    %pixel coordinates: x: ap, + x goes more posterioar (bregma(1)
    %pixel coordinates: y: ml, left brain is at 0, +y goes right, midline
    %at y == 570 (bregma(3)
    pixelCoords = zeros(size(wholeTrack));
    pixelCoords(:,1)=abs(wholeTrack(:,3))/atlasRes+bregma(1); %because were going positive for posterior
    pixelCoords(:,2)=wholeTrack(:,1)/atlasRes + bregma(3);
    pixelCoords(:,3)=wholeTrack(:,2)/atlasRes;
    mo=MEC_Offset/atlasRes;
    
    %av is a matrix containing a regionID for each pixel in 3d volume.
    %Xdim: ap, YDIM: dv, z dim: ml
    %disp(av(round(pixelCoords(1,1)),round(pixelCoords(1,3)),round(pixelCoords(1,2))))

    try
    errorFree = false;
whilecntr = 0;
dAP = [0 0 0];
while ~errorFree && whilecntr<100
    whilecntr = whilecntr +1;
ann = 10;
    start_point = pixelCoords(1,:) + dAP;
    end_point = pixelCoords(2,:)+ dAP;
    p=end_point-start_point;
    p=p/norm(p);
    %p=0.1*vecAllen([3,1,2]); %vec allen is in ml,dv,ap, but pixel coordinats are ap,ml,dv
    %p(3)=abs(p(3));
    try
    isoCtxId = num2str(st.id(strcmp(st.acronym, 'Isocortex')));
    gotToCtx = false;
    cntr = 0;
    while ~(ann==1 && gotToCtx)
        cntr=cntr+1;
        start_point = start_point+p*5; % step 10um, backwards up the track %start point is in ap,ml,dv
        ann = av(round(start_point(1)),round(start_point(3)),round(start_point(2))); %until hitting the top
        %disp(ann)
        if ~isempty(strfind(st.structure_id_path{ann}, isoCtxId))
            % if the track didn't get to cortex yet, keep looking...
            gotToCtx = true;
        end
    end
    errorFree = true;
    catch ME
        dAP=dAP + [-10 0 0];
        errorFree = false;
    end
    
end
    disp(dAP)

    col = 'b';
    catch ME
        %keyboard
        disp(ME.message)
        %disp(av(round(pixelCoords(1,1)),round(pixelCoords(1,2)),round(pixelCoords(1,3))))
        start_point = pixelCoords(2,:);
        plot3(wholeTrack(:,1),wholeTrack(:,3),wholeTrack(:,2),'LineWidth',2,'Color','r')
        col = 'r';
    end
    AllPointsMEC{end+1}=start_point;
    figure(brainfig)
    plot3(pixelCoords(:,1),pixelCoords(:,2),pixelCoords(:,3),'LineWidth',2)
    plot3(start_point(:,1),start_point(:,2),start_point(:,3),'ro')
    figure(controlfig)
    plot3(pixelCoords(:,1),pixelCoords(:,2),pixelCoords(:,3),'LineWidth',1,'Color',col)
    plot3(start_point(:,1),start_point(:,2),start_point(:,3),'ro')
    plot3(pixelCoords(1,1),pixelCoords(1,2),pixelCoords(1,3),'go')
    %plot
    
end
figure(probefig)
xlabel('ML')
ylabel('AP')
zlabel('DV')
axis equal
figure(controlfig)
xlabel('AP')
ylabel('ML')
zlabel('DV')
% figure(brainfig)
% pointList=load('Z:\giocomo\export\data\Projects\AlexA_NP\Histology\AA_190709_1\combined\processed\probe_pointsMEC_L1.mat');
% pointList = pointList.pointList.pointList;
% plot3(pointList{1}(:,3),pointList{1}(:,1),pointList{1}(:,2),'g.')

%%
