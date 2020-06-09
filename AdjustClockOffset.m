function [offset] = AdjustClockOffset(myKsDir)

%myKsDir = '/mnt/analysis/N8/2019-02-07_16-22-02';
filename = fullfile(myKsDir,'messages.events');
fid = fopen(filename,'r');
delimiter = '\t';
startRow = 2;    % ignore first row
textscan(fid, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
offset = cell2mat(textscan(fid,'%f',1));
fclose(fid);

end