function [val, terminate, direction] = event_unstable1(t,z)
    if abs(z(1)) > 0.2 || abs(z(2)) > 0.4
        val = 0;
    else
        val = 1;
    end
    terminate = 1;
    direction = 0;
end