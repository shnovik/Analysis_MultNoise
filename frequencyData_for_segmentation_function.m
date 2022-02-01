function [PopRateSeries] = frequencyData_for_segmentation_function(recordingFilePath, swsFilePath, bindur)
%Converts raw spiking data into a population frequency data, chopping off the first and last bits of each SWS period (less than dt from each side)

% Inputs
    % recordingFilePath -- location of a file from the Levenstein dataset of type recordingname_SSubtypes.mat
    % swsFilePath       -- location of a file from the emailed Levenstein dataset of type recordingname_SlowWaves.events.mat
    % bindur            -- bin duration for binning spikes, in s, typically 0.01

% Outputs
    % PopRateSeries     -- A column vector in which each entry represents a time bin of length bindur, and contains the population frequency
    %                      during that time bin

    % Getting the data:
    load(recordingFilePath)
    load(swsFilePath)
    spike_data = Se_CellFormat; % Se_CellFormat is a variable in the recording file
    decimals = 2;
    sws_times = sws_start_stop_function(swsFilePath);
    first_sws = 1; % If I only want to plot one sws period, first and last can be the same
    last_sws = length(sws_times);

    % Gather all the spike times in one vector:
    all = [];
    for iTrial = 1:length(spike_data)
        spks = cell2mat(spike_data(iTrial)); % Column vector of spiking times for one neuron
        all = [all; spks];
    end


    PopRateSeries = [];
    nbinsL = [];
    test = [];
    for swsRow = first_sws:last_sws
        timePair = sws_times(swsRow,:);
        start = ceil(10^decimals * timePair(1))/10^decimals;
        stop = floor(10^decimals * timePair(2))/10^decimals;
        test = [test; start stop];
        slength = stop - start; % length of current slow wave segment
        nbins = round(slength/bindur); % no actual rounding going on, all my numbers just look like x.0000 and I need them to be x
        nbinsL = [nbinsL; nbins];
    
        part_of_all = [];
        for iSpike = 1:length(all)
            if all(iSpike) >= start & all(iSpike) <= stop
                part_of_all = [part_of_all; all(iSpike)];
            end
        end
    
        h = histogram(part_of_all, nbins);
        part_of_PopRateSeries = (h.Values * (1/length(spike_data)) * (1/bindur)).';
        PopRateSeries = [PopRateSeries; part_of_PopRateSeries];
    
    end

end

