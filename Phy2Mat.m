
% based on example script from N. Steinmetz
% helper function to handle outputs from phy

% defaults
sampleRate = 30000;

%% add the repositories to your path

addpath(genpath('/home/pgupta/Desktop/Code/afterphy'))
addpath(genpath('/home/pgupta/Desktop/Code/spikes'))
addpath(genpath('/home/pgupta/Desktop/Code/npy-matlab'))

%% set paths for where to find your data

myKsDir = '/mnt/data/AZ039/2018-04-23_11-10-05'; %'C:\...\data\myKilosortOutputDirectory';

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

for i = 1:length(sp.cids)
    clear allspikes
    switch sp.cgs(i)
        case 1 % MUA
            
        case 2 % good
            allspikes = sp.st(sp.clu==sp.cids(i)); % all spike times for a given cluster
            
            % get all waveforms for this cluster
            gwfparams.dataDir = myKsDir;                % KiloSort/Phy output folder
            gwfparams.fileName = 'mybinaryfile.dat';    % .dat file containing the raw 
            gwfparams.dataType = 'int16';               % Data type of .dat file (this should be BP filtered)
            gwfparams.nCh = 32;                         % Number of channels that were streamed to disk in .dat file
            gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
            gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out
            gwfparams.spikeTimes =    allspikes*sampleRate; % Vector of cluster spike times (in samples) same length as .spikeClusters
            gwfparams.spikeClusters = sp.cids(i) + 0*gwfparams.spikeTimes; % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
            wf = getWaveForms(gwfparams);
            
            [~,channels] = sort(std(squeeze(wf.waveFormsMean),0,2),'descend');
            channels = ceil(channels/4);
            tetrode = mode(channels(1:4));
            
            
            
            % write spiketimes to .n
            
%             figure;
%             imagesc(squeeze(wf.waveFormsMean))
%             set(gca, 'YDir', 'normal'); xlabel('time (samples)'); ylabel('channel number');
%             colormap(colormap_BlueWhiteRed); caxis([-1 1]*max(abs(caxis()))/2); box off;
            
        otherwise
    end
    
    

end
% plot odor separated raster
%spikeTimes{i,1} = sp.st(sp.clu==sp.cids(i))'; % in seconds


%% Plot spike rasters



% %% Plotting a driftmap
% 
% [spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(myKsDir);
% figure; plotDriftmap(spikeTimes, spikeAmps, spikeDepths);



% %% load synchronization data
% 
% syncChanIndex = 385;
% syncDat = extractSyncChannel(myKsDir, nChansInFile, syncChanIndex);
% 
% eventTimes = spikeGLXdigitalParse(syncDat, lfpFs);
% 
% % - eventTimes{1} contains the sync events from digital channel 1, as three cells:
% % - eventTimes{1}{1} is the times of all events
% % - eventTimes{1}{2} is the times the digital bit went from off to on
% % - eventTimes{1}{2} is the times the digital bit went from on to off
% 
% % To make a timebase conversion, e.g. between two probes:
% % [~,b] = makeCorrection(syncTimesProbe1, syncTimesProbe2, false);
% 
% % and to apply it:
% % correctedSpikeTimes = applyCorrection(spikeTimesProbe2, b);


% %% Looking at PSTHs aligned to some event
% 
% % if you now have a vector of relevant event times, called eventTimes (but
% % not the cell array as above, just a vector):
% 
% window = [-0.3 1]; % look at spike times from 0.3 sec before each event to 1 sec after
% 
% % if your events come in different types, like different orientations of a
% % visual stimulus, then you can provide those values as "trial groups",
% % which will be used to construct a tuning curve. Here we just give a
% % vector of all ones.
% trialGroups = ones(size(eventTimes));
% 
% psthViewer(sp.st, sp.clu, eventTimes, window, trialGroups);
% 
% % use left/right arrows to page through the clusters






