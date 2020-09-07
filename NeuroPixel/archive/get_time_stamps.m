

%%
files = dir('/oak/stanford/groups/giocomo/export/data/Projects/ContrastExperiment_neuropixels/*/*/*/*.nidq.meta');
savedir = '/oak/stanford/groups/giocomo/attialex/timestamps';
imdir = fullfile(savedir,'plots');
if ~isfolder(savedir)
    mkdir(savedir);
    mkdir(imdir);
end
%%
parfor iF=1:numel(files)
    data_dir = files(iF).folder;
    try
    [tsLFP,tsNIDAQ]=get_sync_stamps(data_dir);

    dd=tsLFP-tsNIDAQ;

N = round(numel(dd)/2);

fig = figure('visible','off','Position',[877   510   948   420]);
subplot(1,2,1)
histogram(dd(1:N),'Normalization','probability','BinEdges',[-0.2:0.005:0.2])
hold on
histogram(dd((N+1):end),'Normalization','probability','BinEdges',[-0.2:0.005:0.2])
legend('First Half','second half')
xlabel('difference in timestamp LFP and NIDAQ')
subplot(1,2,2)
plot(tsLFP,dd)
xlabel('recording time in [s]')
ylabel('diff between LFP and NIDAQ')
hold on
x=500;
y=max(dd)/2;
aa=polyfit(tsLFP,dd,1);
text(x,y,sprintf('y = %.4fms *x+%.4fms',aa*1000))
axis tight
sn = strsplit(files(iF).name,'.');
saveas(fig,fullfile(imdir,[sn{1} '.png']))
close(fig)
data_out = matfile(fullfile(savedir,sn{1}),'Writable',true);
data_out.tsLFP = tsLFP;
data_out.tsNIDAQ = tsNIDAQ;
    catch
    end
end
%%
matfiles = dir(fullfile(savedir,'*.mat'));
slopes = nan(numel(matfiles),1);
for iF = 1:numel(matfiles)
    data_out = load(fullfile(matfiles(iF).folder,matfiles(iF).name));
    dd = data_out.tsLFP-data_out.tsNIDAQ;
    aa=polyfit(data_out.tsLFP,dd,1);
    slopes(iF)=aa(1);
end
%%
figure
histogram(slopes*1000)
xlabel('regression slope: [ms/s]')