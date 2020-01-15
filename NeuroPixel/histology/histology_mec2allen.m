histology = process_histology();
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

MEC_Offset = [3.66,2.84,-5.01];
fwireframe = plotBrainGrid([], [], [], true); hold on; 
fwireframe.InvertHardcopy = 'off';
brainfig = gcf;
probefig=figure();
hold on
atlasRes = 0.010; % mm
bregma=allenCCFbregma()
MEC_Offset = [3.0,2.0,-5.25]; %ML,DV,AP %estimate from old paxinos atlas
MEC_angle = 15;
MEC_Offset = [3.3,2.0,-5.0]; % estimated from 2008 allen
for iP= 1:numel(histology.X1)
    mec_entry = [histology.X1(iP),histology.Y1(iP),histology.Z1(iP)]/1000;
    probe_term = [histology.X2(iP),histology.Y2(iP),histology.Z2(iP)]/1000;
    entry_allen = mec_entry([1 3 2]) + MEC_Offset;
    t1=probe_term(1);
    t2=probe_term(3)*cos(deg2rad(15)); %play with this angle
    t3=probe_term(2)*sin(deg2rad(15))+probe_term(2);
    
    term_allen = [t1,t2,t3] + MEC_Offset;
    pos3D=[entry_allen;term_allen];
    vecAllen = -term_allen + entry_allen;
    vecAllen = vecAllen/norm(vecAllen);
    if strcmp(histology.Side{iP},'L')
        %continue
        term_allen(1)=-term_allen(1);
        pos3D(:,1)=-pos3D(:,1);
        vecAllen(1)=-vecAllen(1);
    end
    wholeTrack = [term_allen;term_allen+histology.FinalDepth(iP)*vecAllen/1000];
    figure(probefig)
    scatter3(pos3D(:,1),pos3D(:,3),pos3D(:,2),4,[1 0 0; 0 0 1])
    plot3(wholeTrack(:,1),wholeTrack(:,3),wholeTrack(:,2),'LineWidth',2)
    %plot3(probe_track(:,1),probe_track(:,3),probe_track(:,2))
    pixelCoords = wholeTrack;
    pixelCoords(:,3)=pixelCoords(:,3)*-1;
    pixelCoords = pixelCoords/atlasRes+bregma([3 2 1]);
    
    figure(brainfig)
    plot3(pixelCoords(:,3),pixelCoords(:,1),pixelCoords(:,2),'LineWidth',2)
    %plot3(pixelCoords(1,3),pixelCoords(1,1),pixelCoords(1,2),'ro')
    %plot
    
end
figure(probefig)
xlabel('ML')
ylabel('AP')
zlabel('DV')
axis equal
% figure(brainfig)
% pointList=load('Z:\giocomo\export\data\Projects\AlexA_NP\Histology\AA_190709_1\combined\processed\probe_pointsMEC_L1.mat');
% pointList = pointList.pointList.pointList;
% plot3(pointList{1}(:,3),pointList{1}(:,1),pointList{1}(:,2),'g.')

%%
