function [val, terminate, direction] = event_unstable(t,z)
global start_time
    time_elapsed = clock - start_time;
    if (time_elapsed(5)*60+time_elapsed(6))>20
        val = 0;
    else
        val = 1;
    end
    terminate = 1;
    direction = 0;
end