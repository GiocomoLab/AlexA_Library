matfiles=dir('F:\CatGT2\*\*\imec0_ks2\rez2.mat');

for iF=2:numel(matfiles)
mf = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
figure
plot(mf.rez.dshift)
ylabel('estimated drift')
xlabel('batch #')
saveas(gcf,fullfile(matfiles(iF).folder,'drift_traces.png'))
close(gcf)
end
