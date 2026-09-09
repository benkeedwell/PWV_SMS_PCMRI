function p_out = Holm_Bonferroni(p_in)

%-------------------------------------------------------------------------------%
% APPLY HOLM-BONFERRONI CORRECTION TO P-VALUES
% By Benjamin Keedwell (2026)
%-------------------------------------------------------------------------------%
% INPUTS
% p_in - Uncorrected p-values (N 1)
%-------------------------------------------------------------------------------%
% OUTPUTS
% p_out - Corrected p-values (N 1)
%-------------------------------------------------------------------------------%

% Setup
p_in = p_in(:);
N = length(p_in);
if any(~isfinite(p_in)) || any(p_in < 0 | p_in > 1)
    error('p_in must contain values between 0 and 1');
end

% Sort p-values
[p_sort, idx] = sort(p_in);

% Adjust p-values
p_adj = (N - (1:N)' + 1) .* p_sort;

% Adjusted p-values must remain monotonic
p_adj = cummax(p_adj);

% Adjusted p-values cannot exceed 1
p_adj = min(p_adj, 1);

% Return p-values in original order
p_out = zeros(N, 1);
p_out(idx) = p_adj;

