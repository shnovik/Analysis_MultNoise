
function [emissionSeq] = HMM_format_one_recording(recordingFilePath, swsFilePath, bindur)
% Creates an emmission sequence from spike trains in the input recording to provide as input to an HMM fitting function

% Inputs:
    % recordingFilePath -- location of a file from the Levenstein dataset of type recordingname_SSubtypes.mat
    % swsFilePath       -- location of a file from the emailed Levenstein dataset of type recordingname_SlowWaves.events.mat
    % bindur            -- bin duration used for counting the number of spikes per bin, measured in seconds (typically bindur = 0.01 s)

% Outputs:
    % emissionSeq       -- cell array in which each cell represents one period of slow wave sleep. Each cell contains a matrix with one row
    %                      per neuron/channel and one column per time bin. Each entry is number of times that neuron spiked in that time bin.

    
    load(recordingFilePath)
    load(swsFilePath)
        
    spike_data = Se_CellFormat; % Se_CellFormat is a cell array variable found in recordingname_SSubtypes.mat
    decimals = 2; % this is used in rounding starts/stops correctly so that there's an integer number of bins
    sws = sws_start_stop_function(swsFilePath);
    first_sws = 1;
    last_sws = length(sws);
    
    emissionSeq = cell(1, length(sws));

    
    for swsRow = first_sws:last_sws
        timePair = sws(swsRow,:);
        start = ceil(10^decimals * timePair(1))/10^decimals;
        stop = floor(10^decimals * timePair(2))/10^decimals;
        slength = stop - start;
        nbins = round(slength/bindur); % no actual rounding going on, all my numbers just look like x.0000 and I want them to be x
        swsPeriodMatrix = zeros(length(spike_data), nbins);
        
        for iTrial = 1:length(spike_data)
            spike_train = cell2mat(spike_data(iTrial)); % Column vector of spiking times for one neuron
            
            sws_spks = [];   
            for iSpike = 1:length(spike_train)
                if spike_train(iSpike) >= start & spike_train(iSpike) <= stop
                    sws_spks = [sws_spks; spike_train(iSpike)];
                end
            end
            
            histogramValues = histcounts(sws_spks, nbins);
            swsPeriodMatrix(iTrial,:) = histogramValues;
            
            %h = histogram(sws_spks, nbins);
            %swsPeriodMatrix(iTrial,:) = h.Values;
                 
        end
    
        emissionSeq{swsRow} = swsPeriodMatrix;
        
    end
