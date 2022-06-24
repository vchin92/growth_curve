% Gibbs sampling for auxiliary variable for random truncation
function [output] = gibbs_u(grp, weight)

output = rand(length(grp), 1) .* weight(grp);

end