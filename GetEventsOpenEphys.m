% script to get events from the openephys session
addpath(genpath('/home/pgupta/Desktop/Code/open-ephys-analysis-tools'))
addpath(genpath('/home/pgupta/Desktop/Code/mouse-lever-task'))

datadir = '/mnt/data/N8/2019-01-26_19-24-28/';
filename = fullfile(datadir,'all_channels.events');
behavior_filename = fullfile(datadir,'N8_20190126_r0.mat');

[events, timestamps, info] = load_open_ephys_data(filename);
% trial ON-OFF events are ID = 3
TrialOn = timestamps(intersect(find(info.eventId),find(info.eventType==3)));
TrialOff = timestamps(intersect(find(~info.eventId),find(info.eventType==3)));

[TrialInfo] = GetBehaviorEvents(behavior_filename);

