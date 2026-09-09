function x = BeeSwarm(y, x0, ry, rx)
%-------------------------------------------------------------------------------%
% JITTER X-COORDS FOR BEESWARM PLOT
% By Benjamin Keedwell (2026)
%-------------------------------------------------------------------------------%
% Introduce jitter for beeswarm x-coordinates.
%-------------------------------------------------------------------------------%
% INPUTS
% y - y-coordinates (N 1)
% x0 - x-coordinate (Scalar)
% ry - Minimum y-spacing
% rx - Minimum x-spacing
%-------------------------------------------------------------------------------%
% OUTPUTS
% x - Jittered x-coordinates (N 1)
%-------------------------------------------------------------------------------%

% Setup
[ys, idx] = sort(y);

N = length(y);

x = zeros(N, 1);

% Iterate through sorted datapoints
for i = 1:N
    xi = 0;
    k = 0;

    while true
        overlap = false;

        % Compare with previous points
        for j = 1:i-1
            dy = abs(ys(i) - ys(j));
            if dy < ry
                dx = abs(xi - x(j));
                if dx < rx
                    overlap = true;
                    break
                end
            end
        end

        % For no overlap, stop shifting position
        if ~overlap
            break
        end

        % Alternate between left/right shifts
        % (i+k) stops datapoints from mostly being shifted negatively
        k = k + 1;
        xi = (-1)^(i+k) * ceil(k/2) * rx;

        if k == 100
            error('Beeswarm jitter has ran 100 iterations')
        end

    end

    % Final x-position
    x(i) = xi;
end

% Restore original order
xtemp = x;
x(idx) = xtemp;

% Shift to input x-position
x = x + x0;

end