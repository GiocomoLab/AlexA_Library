% GLM_run_model_on_all_cells.m
% Script to run Kiah's GLM on neuropixels data.
% MGC 2/12/2019

%% data and params

% make sure paths are correct
restoredefaultpath
addpath(genpath('C:\Users\giocomolab\Dropbox\Work\neuropixels\functions'));
addpath(genpath('C:\Users\giocomolab\Dropbox\Work\neuropixels\spikes'));
addpath(genpath('C:\Users\giocomolab\Dropbox\Work\neuropixels\spline-lnp-model-VR'));

% where to find data
data_dir = 'C:\Users\giocomolab\Dropbox\Work\neuropixels\data\';
session_name = 'npF3_1019_contrasttrack_gainchanges_contrast_1';

% some params
params = readtable('UniversalParams.xlsx');

%% load data

load(fullfile(data_dir,strcat(session_name,'.mat')));
speed = calcSpeed(posx,params);
% map negative speeds to zero
speed(speed<0) = 0;
% take position mod length of track (AFTER computing speed)
posx = mod(posx,params.TrackEnd); 
good_cells = sp.cids(sp.cgs==2);
numcells = numel(good_cells);

%% Run glm on all cells;

bestModels = cell(numcells,1);
allModelTestFits = cell(numcells,1);
final_pval = nan(numcells,1);
for i = find(good_cells==921)
    
    fprintf('\nCell %d/%d\n',i,numcells);
    
    % spike times
    spike_t = sp.st(sp.clu==good_cells(i));
    % spiketrain
    spiketrain = histc(spike_t,post);
    [~,~,spike_idx] = histcounts(spike_t,post);

    % run glm
    [bestModels{i}, allModelTestFits{i}, ~, final_pval(i), fig1] = create_glm(posx,speed,spiketrain,params);
    % plot spike raster
    subplot(4,3,[2 5 8]);
    plot(posx(spike_idx),trial(spike_idx),'k.');
    ylim([0 max(trial)+1]); xlim([params.TrackStart params.TrackEnd]);
    xlabel('pos'); ylabel('trial');
    % save fig and close
    %saveas(fig1, sprintf('images/glm/%s/%d.png',session_name,good_cells(i)), 'png'); 
    %close(1);
end