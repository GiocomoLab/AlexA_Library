bin_path = 'F:\AA\output\20210622_M217\vr_baseline-004\output\suite2p\plane1\data.bin'
finfo=dir(bin_path);
nframes = finfo.bytes/(512*512*2);
fi_mem = memmapfile(bin_path,'Format',{'int16',[512 512 nframes],'x'});
%%
view_stack(fi_mem.Data.x,'map','jet')