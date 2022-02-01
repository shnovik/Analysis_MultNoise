function Raster_Frequency_Plots(recordingFilePath, swsFilePath, first_sws, last_sws, down_pairs, up_pairs, bindur)
% Creates a raster plot and a population frequency histogram, with Up and Down states shaded in different colors

% Inputs
    % recordingFilePath -- file path to a file from the Levenstein dataset of type recordingname_SSubtypes.mat
    % swsFilePath       -- file path to a file from the emailed Levenstein dataset of type recordingname_SlowWaves.events.mat
    % first_sws         -- the index of the first period of slow wave sleep to plot, e.g. '1' to plot the first period
    % last_sws          -- the index of the last period of slow wave sleep to plot (same as first_sws if only plotting one period of slow wave
                            % sleep)
    % down_pairs        -- start and stop times of down states (something like either recordingname_SlowWaves.ints.DOWN from Levenstein data or
                            % downTimes from HMM_stateTimestamps.m output
    % up_pairs          -- start and stop times of up states
    % bindur            -- duration of bins for the frequency histogram in s, typically 0.01

    % Initializing the figure:
    figure()
    ax = subplot(4, 1, 1);

    % Getting the data:
    load(recordingFilePath)
    load(swsFilePath)
    data = Se_CellFormat; % a variable in the recording file
    sws_times = sws_start_stop_function(swsFilePath);
    xstart = sws_times(first_sws, 1);
    xstop = sws_times(last_sws, 2);
    all = []; % this is used in making the frequency histogram


    % -------------
    % Raster plot
    % -------------
    
    ax.XLim = [xstart xstop];
    %ax.XLim = [xstart xstart+15];
    ax.YLim = [0 length(data)];
    ax.XLabel.String = "Time [s]";
    ax.YLabel.String = "Neurons";

    hold on

    % shading raster plot based on HMMs and Levenstein segmentations
    for down = 1:length(down_pairs)
        patch([down_pairs(down, 1) down_pairs(down, 2) down_pairs(down, 2) down_pairs(down, 1)], [max(ylim) max(ylim) 0 0], [0.3010, 0.7450, 0.9330])
    end
    hold on
    for up = 1:length(up_pairs)
        patch([up_pairs(up, 1) up_pairs(up, 2) up_pairs(up, 2) up_pairs(up, 1)], [max(ylim) max(ylim) 0 0], [0.9290, 0.6940, 0.1250])
    end
    hold on

    % plotting the spikes
    for iTrial = 1:length(data)
    
        spks = cell2mat(data(iTrial)).'; % Row vector of spiking times for one neuron
    
        % Get only the spikes occuring during the slow wave sleep periods of interest:
        sws_spks = [];
        for iSpike = 1:length(spks)
            for swsRow = first_sws:last_sws
                timePair = sws_times(swsRow,:);
                if spks(iSpike) > timePair(1) & spks(iSpike) < timePair(2)
                    sws_spks = [sws_spks spks(iSpike)];
                end
            end
        end
    
        xspikes = repmat(sws_spks, 2, 1);       % spks repeated 2 times, for plotting      
        yspikes = nan(size(xspikes));           % to be filled with start and end y values for each mark on the plot
    
        if ~isempty(yspikes)
            yspikes(1,:) = iTrial - 1;
            yspikes(2,:) = iTrial;
        end

        plot(xspikes, yspikes, 'Color', 'k')
    
        all = horzcat(all, sws_spks); % Concatenate all spikes from all neurons into all, for frequency histogram
    end
    
    
    
    % -------------
    % Frequency plot
    % -------------
    
    slength = xstop - xstart; % length of slow wave segment being plotted

    ax = subplot(4, 1, 2);
    ax.XLim = [xstart xstop];
    %ax.XLim = [xstart xstart+15];
    ax.YLim = [0 5];


    % Shading the frequency histogram based on Up/Down segmentation

    % code used to do this for Ann Christin's segmentation
    %crossings = [ 122,  167,  171,  178,  181,  310,  312,  606,  609, 1073, 1080, 1441, 1442, 1754, 1756, 1868, 1879, 1886, 1892, 2216, 2218]; % Get this from running mov_avg_crossings_durations_2thresholds in changedStatistics.py
    %crossing_times = zeros(0, length(crossings));
    %for i = 1:length(crossings)
    %    crossing_times(i) = (crossings(i) * bindur) + xstart; %sws(first_sws, 1); % + xstart; 
    %end

    %if rem(length(crossing_times), 2) == 0
    %    num_downs = (length(crossing_times))/2;
    %    num_ups = num_downs - 1;
    %else
    %    num_downs = (length(crossing_times) - 1)/2; 
    %    num_ups = num_downs;
    %end

    %down_pairs = zeros(num_downs, 2);
    %up_pairs = zeros(num_downs - 1, 2);
    %for i = 1:num_downs
    %    downCross = 1 + (i-1)*2;
    %    down_pairs(i, 1) = crossing_times(downCross);
    %    down_pairs(i, 2) = crossing_times(downCross+1);
    %end
    %for i = 1:num_ups
    %    upCross = i*2;
    %    up_pairs(i, 1) = crossing_times(upCross);
    %    up_pairs(i, 2) = crossing_times(upCross + 1);
    %end

    % shading based on HMMs and Levenstein segmentations
    for down = 1:length(down_pairs)
        patch([down_pairs(down, 1) down_pairs(down, 2) down_pairs(down, 2) down_pairs(down, 1)], [max(ylim) max(ylim) 0 0], [0.3010, 0.7450, 0.9330])
    end
    hold on
    for up = 1:length(up_pairs)
        patch([up_pairs(up, 1) up_pairs(up, 2) up_pairs(up, 2) up_pairs(up, 1)], [max(ylim) max(ylim) 0 0], [0.9290, 0.6940, 0.1250])
    end

    
    % plot the frequency histogram
    nbins = round(slength/bindur);
    h = histogram(all, nbins);
    h.FaceColor = 'k';
    %avFreq = h.Values * (1/length(data)) * (1/bindur); % for use in Ann Christin's python segmentation 

    ax.XLabel.String = "Time [s]";
    ax.YLabel.String = 'Spikes/Bin';
    %ax.XLim = [xstart xstop];

    % Convert y axis from spikes/bin to spikes/s
    newlabel = [];
    for iLab = 1:length(ax.YTickLabel)
        lab = str2num(ax.YTickLabel{iLab});
        conv = (lab / length(data)) * (1/bindur);
        newlabel = [newlabel, roundn(conv, -1)];
    end

    ax.YTickLabel = newlabel;
    ax.YLabel.String = 'Firing rate [Hz]';

end

