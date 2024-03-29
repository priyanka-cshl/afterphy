function [offset] = AdjustClockOffset(myKsDir)

% %myKsDir = '/mnt/analysis/N8/2019-02-07_16-22-02';
% filename = fullfile(myKsDir,'messages.events');
% fid = fopen(filename,'r');
% delimiter = '\t';
% startRow = 2;    % ignore first row
% textscan(fid, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
% offset = cell2mat(textscan(fid,'%f',1));
% fclose(fid);
% 
% if isempty(offset) % offset is in the first line
%     fid = fopen(filename,'r');
%     delimiter = '\t';
%     startRow = 1;    % ignore first row
%     textscan(fid, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
%     offset = cell2mat(textscan(fid,'%f',1));
% end
% 
% fclose(fid);

if ~isempty(dir(fullfile(myKsDir,'*_ADC1.continuous')))
    foo = dir(fullfile(myKsDir,'*_ADC1.continuous')); % pressure sensor
else
    numChannels = numel(dir(fullfile(myKsDir,'1*_*.continuous')));
    foo = dir(fullfile(myKsDir,['1*_', num2str(numChannels-7) ,'.continuous'])); % pressure sensor
end
filename = fullfile(myKsDir,foo.name);
[Auxdata1, timestamps, ~] = load_open_ephys_data(filename); % data has channel IDs
offset = timestamps(1);

end