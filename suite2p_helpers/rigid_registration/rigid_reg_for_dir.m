function rigid_reg_for_dir(input_dir)

sbx_files = dir(fullfile(input_dir,'*.sbx'));

for ii=1:length(sbx_files)
    fprintf('Now working on: %s \n',sbx_files(ii).name)
    [~,fn,~]=fileparts(sbx_files(ii).name);
    sbxread(fullfile(input_dir,fn),1,1);            % read one frame to read the header of the image sequence
    global info;
    fprintf('Loading data \n')
    data=sbxread(fullfile(input_dir,fn),0,info.max_idx);
    %data=sbxread(fullfile(input_dir,fn),0,2000);
    for iL=1:2
        fprintf('Running registration for layer %d \n',iL)
        [dx,dy,template]=rigid_registration(squeeze(data(1,:,99:700,iL:2:end)));
        dX{iL}=dx;
        dY{iL}=dy;
        templates{iL}=template;
        fig=figure;
        imagesc(template)
        saveas(fig,fullfile(input_dir,sprintf([fn 'reg_L%d.png'],iL )))
        close(fig)
    end
    registration.dX=dX;
    registration.dY=dY;
    registration.templates=templates;

    save(fullfile(input_dir,fn),'registration','-append')


end
