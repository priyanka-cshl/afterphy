
% based on example script from N. Steinmetz
% helper function to handle outputs from phy

%% add the repositories to your path

addpath(genpath('/home/pgupta/Desktop/Code/afterphy'))
addpath(genpath('/home/pgupta/Desktop/Code/spikes'))
addpath(genpath('/home/pgupta/Desktop/Code/npy-matlab'))
addpath(genpath('/home/pgupta/Desktop/Code/matlab_utils/plotSpikeRaster'))
addpath(genpath('/home/pgupta/Desktop/Code/open-ephys-analysis-tools'))
addpath(genpath('/home/pgupta/Desktop/Code/mouse-lever-task'))

%% set paths for where to find your data

myKsDir = '/mnt/data/N8/2019-01-26_19-24-28'; %'C:\...\data\myKilosortOutputDirectory';

%myEventTimes = load('C:\...\data\someEventTimes.mat'); % a vector of times in seconds of some event to align to
% load events from the open ephys data
filename = fullfile(myKsDir,'all_channels.events');
[events, timestamps, info] = load_open_ephys_data(filename);
% trial ON-OFF events are ID = 3
TrialOn = timestamps(intersect(find(info.eventId),find(info.eventType==3)));
TrialOff = timestamps(intersect(find(~info.eventId),find(info.eventType==3)));

% get odor IDs from behavior file
behavior_dir = '/mnt/data/N8/2019-01-26_19-24-28/';
behavior_filename = fullfile(behavior_dir,'N8_20190126_r0.mat');
[TrialInfo, Traces] = GetBehaviorEvents(behavior_filename);

%% Loading data from kilosort/phy easily

sp = loadKSdir(myKsDir);

% sp.st are spike times in seconds (for all spikes)
% sp.clu are cluster identities (for all spikes)
% sp.cids is list of unqiue clusters
% sp.cgs are cluster defs (1 = MUA, 2 = good, 3 = Unsorted??) (1/cluster)
% spikes from clusters labeled "noise" have already been omitted

%% Split data by clusters and by trials
recordinglength = max(sp.st);
clear rasters

Lever = Traces.Lever;
% align to trial off
for t = 1:size(Lever,1)
    MyLever(t,:) = circshift(Lever(t,:),TrialInfo.Duration(t));
end

for i = 12:length(sp.cids)
    switch sp.cgs(i)
        case 2
            
        otherwise
    end
    allspikes = sp.st(sp.clu==sp.cids(i));
    clear spikeTimes
    clear validspikeTimes
    % parse to trials
    for j = 1:min(numel(TrialOn),numel(TrialOff))
        thisTrialSpikes = allspikes(find(allspikes>=TrialOn(j) & allspikes<=TrialOff(j)));
        
        % do this again - delete any spikes that occur when not in TZ
        validspikes = [];
        for foo = 1:numel(TrialInfo.StayTimeStart{j})
            zstart = TrialOn(j) + 0.002*TrialInfo.StayTimeStart{j}(foo);
            zstop = TrialOn(j) + 0.002*(TrialInfo.StayTimeStart{j}(foo) + TrialInfo.StayTime{j}(foo));
            validspikes = [validspikes; find(thisTrialSpikes>=zstart & thisTrialSpikes<=zstop)];
        end
        
        thisTrialSpikes = thisTrialSpikes - TrialOff(j);
        spikeTimes{j,1} = thisTrialSpikes';
        if size(thisTrialSpikes(validspikes)',1) == 0
            validspikeTimes{j,1} = validspikeTimes{4,1};
        else
            validspikeTimes{j,1} = thisTrialSpikes(validspikes)';
        end
        
        
    end
    for k = 1:3 % 3 odors
        trial_sequence = [];
        trial_marks = [];
        for m = 1:12
            f = intersect(find(TrialInfo.Odor==k),find(TrialInfo.TargetZoneType==m));
            trial_sequence = [trial_sequence; f];
            trial_marks = [trial_marks length(trial_sequence)];
        end
        raster{i,k} = spikeTimes(trial_sequence);
        
        subplot(3,1,k);
        LineFormat.Color = 'b';
        plotSpikeRaster(raster{i,k},'PlotType','vertline','LineFormat',LineFormat);
        
        LineFormat.Color = 'r';
        plotSpikeRaster(validspikeTimes(trial_sequence),'PlotType','vertline','LineFormat',LineFormat);
        
        for foo = 1:numel(trial_marks)
            line(get(gca,'XLim'), trial_marks(foo)*[1 1], 'color', 'k')
        end
        
    end
end
% plot odor separated raster
%spikeTimes{i,1} = sp.st(sp.clu==sp.cids(i))'; % in seconds


%% Plot spike rasters



%% Plotting a driftmap

[spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(myKsDir);
figure; plotDriftmap(spikeTimes, spikeAmps, spikeDepths);



%% load synchronization data

syncChanIndex = 385;
syncDat = extractSyncChannel(myKsDir, nChansInFile, syncChanIndex);

eventTimes = spikeGLXdigitalParse(syncDat, lfpFs);

% - eventTimes{1} contains the sync events from digital channel 1, as three cells:
% - eventTimes{1}{1} is the times of all events
% - eventTimes{1}{2} is the times the digital bit went from off to on
% - eventTimes{1}{2} is the times the digital bit went from on to off

% To make a timebase conversion, e.g. between two probes:
% [~,b] = makeCorrection(syncTimesProbe1, syncTimesProbe2, false);

% and to apply it:
% correctedSpikeTimes = applyCorrection(spikeTimesProbe2, b);


%% Looking at PSTHs aligned to some event

% if you now have a vector of relevant event times, called eventTimes (but
% not the cell array as above, just a vector):

window = [-0.3 1]; % look at spike times from 0.3 sec before each event to 1 sec after

% if your events come in different types, like different orientations of a
% visual stimulus, then you can provide those values as "trial groups",
% which will be used to construct a tuning curve. Here we just give a
% vector of all ones.
trialGroups = ones(size(eventTimes));

psthViewer(sp.st, sp.clu, eventTimes, window, trialGroups);

% use left/right arrows to page through the clusters


%% PSTHs across depth

depthBinSize = 80; % in units of the channel coordinates, in this case �m
timeBinSize = 0.01; % seconds
bslWin = [-0.2 -0.05]; % window in which to compute "baseline" rates for normalization
psthType = 'norm'; % show the normalized version
eventName = 'stimulus onset'; % for figure labeling

[timeBins, depthBins, allP, normVals] = psthByDepth(spikeTimes, spikeDepths, ...
    depthBinSize, timeBinSize, eventTimes, window, bslWin);

figure;
plotPSTHbyDepth(timeBins, depthBins, allP, eventName, psthType);


%% Loading raw waveforms

% To get the true waveforms of the spikes (not just kilosort's template
% shapes), use the getWaveForms function:

gwfparams.dataDir = myKsDir;    % KiloSort/Phy output folder
apD = dir(fullfile(myKsDir, '*ap*.bin')); % AP band file from spikeGLX specifically
gwfparams.fileName = apD(1).name;         % .dat file containing the raw
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 385;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out
gwfparams.spikeTimes = ceil(sp.st(sp.clu==155)*30000); % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = sp.clu(sp.clu==155);

wf = getWaveForms(gwfparams);

figure;
imagesc(squeeze(wf.waveFormsMean))
set(gca, 'YDir', 'normal'); xlabel('time (samples)'); ylabel('channel number');
colormap(colormap_BlueWhiteRed); caxis([-1 1]*max(abs(caxis()))/2); box off;
