% Gibbs sampling for concentration parameter
function [output] = gibbs_eta(grp, eta, prior_eta)

v     = prior_eta.v;
delta = prior_eta.delta;

x      = betarnd(eta, length(grp));
output = gamrnd(v + length(unique(grp)), 1 / (delta - log(x)));

end