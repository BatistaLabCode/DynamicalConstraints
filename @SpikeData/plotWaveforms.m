function [vpp,F] = plotWaveforms(S,varargin)
% Plot waveforms method for SpikeData object
%
% This function plots the waveforms in the provided SpikeData object.

% Create plot and axes handles
F = figure('Position',[10 10 1000 1000]);
H = genArrayAxes();

% Optional arguments
sortCode = 1:5;         % Valid sort codes
dsFrac = .01;          % Fraction of waveforms to plot for each channel/sort code
YLim = 'auto';          % Specify if scale should be automatic

assignopts(who,varargin);

% Get channel and sort information
channel = [S.channel];
sort = [S.sort];
uniCh = unique(channel);
nCh = length(uniCh);
vpp = nan(nCh,1);

% Loop over channels
for i = 1:nCh
    
    % Find all waveforms with the appropriate sort code(s)
    chMask = (channel == uniCh(i));
    srtMask = ismember(sort,sortCode);
    mask = chMask & srtMask;
    w = [S(mask).waveform];
   
    % Calculate peak-to-peak voltage for all waveforms.  Instead of
    % reporting the peak-to-peak value of the average, instead report the
    % 95% upper bound -- this should capture the largest units on the
    % electrode while still removing any outliers.
    nSpikes = size(w,2);
    if nSpikes > 0
        vpp_temp = max(w) - min(w);
        vpp(i) = quantile(vpp_temp, 0.95) * 1e6; % Convert to uV;
    end
    
    % Downsample (plotting purposes only)
    wAvg = mean(w,2); % Calcualte average using all spikes
    spikeInd = randperm(nSpikes);
    nDs = round(nSpikes * dsFrac);
    spikeInd = spikeInd(1:nDs);
    w = w(:,spikeInd);
    
    % Plot
    axes(H(i))
    hold on
    for j = 1:nDs
        plot(w(:,j),'color',ones(1,3)*.65)
    end

    if nDs > 0
        % Plot average waveform and peak-to-peak voltage
        plot(wAvg,'k','LineWidth',2)
        t = sprintf(' = %0.1f ',vpp(i));

        % Scale y-axis limits
        axis tight
        
        if strcmp(YLim,'auto')
            yLim = get(gca,'YLim');
            dY = yLim(2) - yLim(1);
            yLim(1) = yLim(1) - dY*.1;
            yLim(2) = yLim(2) + dY*.1;
        else
            yLim = YLim / 1e6;  % Convert from uV
        end
        
        set(H(i),'xLim',[1 size(w,1)],'yLim',yLim)
        yLim = get(gca,'YLim');
        text(2,yLim(2),['V_{pp}' t '{\mu}V'],'color','r', ...
            'VerticalAlignment','top','FontWeight','bold')
    end
    
    set(H(i),'XTickLabel',[],'YTickLabel',[])
    drawnow
end