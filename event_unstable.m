function [val, terminate, direction] = event_unstable(t,z)
    if abs(z(1)) > 10
        val = 0;
    else
        val = 1;
    end
    terminate = 1;
    direction = 0;
end