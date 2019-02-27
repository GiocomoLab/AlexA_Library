function rigid_reg_for_dir(input_dir)

sbx_files = dir(fullfile(input_dir,'*.sbx'));

for ii=1:length(sbx_files)
    fprintf('Now working on: %s \n',sbx_files(ii).name)
    [~,fn,~]=fileparts(sbx_files(ii).name);
    sbxread(fullfile(input_dir,fn),1,1)            % read one frame to read the header of the image sequence
    global info;
    data=sbxread(fullfile(input_dir,fn),0,info.max_idx);
    for iL=1:2
        [dx,dy,template]=rigid_registration(squeeze(data(1,:,100:690,iL:2:end)));
        dX{iL}=dx;
        dY{iL}=dy;
        templates{iL}=template;
    end
    registration.dX=dX;
    registration.dY=dY;
    registration.templates=templates;
    save(fullfile(input_dir,fn),registration,'-append')


end
