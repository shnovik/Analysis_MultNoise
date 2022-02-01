
function [sws_times] = sws_start_stop_function(slow_waves_file)
% Generates a list of start and stop times for overall slow wave sleep periods

% Inputs
    % slow_waves_file -- file path to a file from the emailed Levenstein dataset of type recordingname_SlowWaves.events.mat
    
% Outputs
    % 2*n matrix containing start and stop times of slow wave periods, n = total number of slow wave periods

% Get data
load(slow_waves_file);
up_data = SlowWaves.ints.UP;
down_data = SlowWaves.ints.DOWN;

% Initialize slow wave timestamp array
sws_times = [];

% Set first slow wave start
[start_sws, current_state] = start_helper(1, 1, up_data, down_data);

% Set termination point
if length(up_data) < length(down_data)
    terminate = length(up_data);
else
    terminate = length(down_data);
end



% Find all the slow wave start and stops

u = 1; %row of up data
d = 1; %row of down data
while u < terminate
    if current_state == "up"
        if down_data(d, 1) == up_data(u, 2)
            current_state = "down";
        else
            stop_sws = up_data(u, 2);
            sws_times = [sws_times ; start_sws stop_sws];
            [start_sws, current_state] = start_helper(u+1, d, up_data, down_data);
        end
        u = u + 1;
    else
        if up_data(u, 1) == down_data(d, 2)
            current_state = "up";
        else
            stop_sws = down_data(d, 2);
            sws_times = [sws_times ; start_sws stop_sws];
            [start_sws, current_state] = start_helper(u, d+1, up_data, down_data);
        end
        d = d+1;
    end

end

