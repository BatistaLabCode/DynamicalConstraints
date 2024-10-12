% combineSorts      Combine sort code information method for SpikeData
% class.
%
% Usage:
%   [S] = S.combineSorts;
%
% Author:   Alan D. Degenhart
% Created:  2017.12.19
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function [S] = combineSorts(S,varargin)
% Optional arguments
sortCodes = [1:30];     % Sort codes to merge

assignopts(who,varargin);

% Get channel and sort information
chn = [S.channel];
srt = [S.sort];
uniCh = unique(chn);
nUniCh = length(uniCh);
removeMask = zeros(1,length(S));

% Loop over unique channels
Scomb = repmat(SpikeData(),1,nUniCh);
for i = 1:nUniCh
    % Get spike data to merge for given channel
    chnMask = ismember(chn,uniCh(i));
    srtMask = ismember(srt,sortCodes);
    mask = chnMask & srtMask;
    Stmp = S(mask);
    
    % Update mask of spike data elements to remove
    removeMask = removeMask | mask;
    
    % Get all spike times and waveforms and sort
    st = [Stmp.spikeTimes];
    wvf = [Stmp.waveform];
    [st,I] = sort(st,'ascend');
    
    % Sort waveforms if available
    if S(1).includesWaveforms
        wvf = wvf(:,I);
    end
    
    % Add sort codes and waveform to new SpikeData object
    Scomb(i).channel = uniCh(i);
    Scomb(i).sort = sortCodes(1);
    Scomb(i).source = S(1).source;
    Scomb(i).timeStampUnits = S(1).timeStampUnits;
    Scomb(i).spikeTimes = st;
    Scomb(i).waveform = wvf;
    Scomb(i).includesWaveforms = S(1).includesWaveforms;
end

% Remove spike data objects that have been combined and add combined
% objects
S = S(~removeMask);
S = [S Scomb];

% Sort spike data objects by channel and sort code
chn = [S.channel];
Ssrt = repmat(SpikeData(),1,length(S));
onset = 1;
for i = 1:nUniCh
    % Get all sorts for current channel
    chMask = ismember(chn,uniCh(i));
    Stmp = S(chMask);
    
    % Get sort code and sort
    srt = [Stmp.sort];
    [~,I] = sort(srt,'ascend');
    
    % Add data to sorted spike data structure
    offset = onset + length(srt) - 1;
    Ssrt(onset:offset) = Stmp(I);
    onset = offset + 1;
end 

% Update output sorted SpikeData array
S = Ssrt;