function proj_meta = get_proj_meta(exp_table,root)
site_ctr = 1;
for iR = 1:height(exp_table)
    try
        animal=exp_table.Animal{iR};
        site = exp_table.Site(iR);
        exp = exp_table.Exp(iR);
        exp_string_sbx = sprintf([animal '_%03d_%03d'],site,exp);
        exp_string_other = sprintf([animal '_%d_%d'],site,exp);
        
        
        
        animal_root = fullfile(root,animal);
        
        sbxfile = fullfile(animal_root,exp_string_sbx);
        
        sbxread(sbxfile,0,1);
        global info
        filename_ROI = fullfile(sbxfile,'suite2p','combined');
        dF=calcdF_s2p(filename_ROI);
        
        logfiles = dir([animal_root , '\*.log']);
        for iF=1:length(logfiles)
            if contains(logfiles(iF).name,exp_string_other)
                filename_log=logfiles(iF).name;
            end
        end
        
        vr_data = sync_2p_vr(fullfile(animal_root,filename_log),info);
        clear info
        proj_meta(site_ctr).animal = animal;
        proj_meta(site_ctr).site = site;
        proj_meta(site_ctr).exp = exp;
        proj_meta(site_ctr).dFF = dF;
        proj_meta(site_ctr).vr_data=vr_data;
        
        site_ctr = site_ctr+1;
        
    catch ME
        disp(ME.identifier)
        fprintf('Failed for Animal, %s, row number %d \n',animal,iR)
    end
end

end
