function [start_sws, current_state] = start_helper(u, d, up_data, down_data)
%{
Helper function for sws_start_stop.m

After the end of an sws period is found, this function is used to
determine whether the next sws period starts with the next up state or
the next down state. It sets the timestamp at which the next sws period
starts and determines which state this is.

Inputs:
   up_data and down_data: Data of start and stop times for up and down states
   u and d: the indeces of the next possible up state and down state
%}

if up_data(u, 1) < down_data(d, 1)
    start_sws = up_data(u, 1);
    current_state = "up";
else
    start_sws = down_data(d, 1);
    current_state = "down";
end
end

