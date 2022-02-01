function [upTimes, downTimes] = HMM_stateTimestamps(emissionSeq, estTrans, estEmis, estPi0, swsFile, binsize)
% Given an emission sequence, an HMM, and a file from which sws period
% start and stop times can be determiend, finds the start and stop times (in s) of all Down
% states and of all Up states for the recording

% Inputs
    % emissionSeq -- An output of HMM_format_one_recording.m
    % estTrans    -- Transition matrix between states, an output of fitHMM.m
    % estEmis     -- Emission matrix, an output of fitHMM.m
    % estPi0      -- Initial probability distribution, an output of fitHMM.m
    % swsFile     -- Location of a file from the emailed Levenstein dataset of type recordingname_SlowWaves.events.mat
    % binsize     -- Bins used to create emissionSeq, typically 0.01 s
    
% Outputs
    % upTimes     -- Start and stop times of every Up state, in s
    % downTimes   -- Start and stop times of every Down state, in s

upTimes = [];
downTimes = [];

decimals = 2;

numTrial = length(emissionSeq);

for iTrial = 1:numTrial
    states{iTrial} = hmmviterbiPoisson(emissionSeq{iTrial},estTrans,estEmis, estPi0); %get the state sequence for each period of slow wave sleep
end;

swsTimes = sws_start_stop_function(swsFile);

% Convert the Viterbi state sequence into start/stop timestamps
for iTrial = 1:numTrial
    start = ceil(10^decimals * swsTimes(iTrial, 1))/10^decimals;
    timeDif = 0;
    one_stateSeq = states{iTrial};
    for state_part = 1:length(one_stateSeq) - 1
        timeDif = timeDif + binsize;
        if one_stateSeq(state_part) ~= one_stateSeq(state_part + 1) || state_part == length(one_stateSeq) - 1
            stop = start + timeDif;
            if one_stateSeq(state_part) == 1
                upTimes = [upTimes; start stop];
            else
                downTimes = [downTimes; start stop];
            end
            start = stop;
            timeDif = 0;
        end
    end

end

