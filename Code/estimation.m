%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Article : Multiclass classification of growth curves using random change 
%            points and heterogeneous random effects
%
%  Authors : Vincent Chin, Jarod Lee, Louise M. Ryan, Robert Kohn, and
%            Scott A. Sisson
%
%  Date : 11.06.2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data: df_fixed.mat/df_random.mat
load('df_random.mat');

% Replication number (1-100)
nrep = 1;

% Number of individuals modulo 50 (maximum value of 8)
% 1 - 50 individuals, 2 - 100 individuals etc.
nsize = 3;

% Indicator for variable knot
% 0: fixed knot
% 1: variable knot
variable_knot = 1;

% Number of knots
nknot = 2;

% Number of MCMC iterations
n_mcmc = 100000;

% Number of burnin
burnin = n_mcmc/2;

% Thinning
thin   = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prior distribution parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_pred_xf      = 1;
prior_alpha.mat  = 25 * eye(num_pred_xf);
prior_alpha.mean = zeros(num_pred_xf, 1);

variance_scale  = 5;

prior_var.delta = 1;
prior_var.A     = variance_scale;
prior_var.a     = 1;
prior_var.Xi    = 1;

prior_rint.delta = 1;
prior_rint.scale = variance_scale;

num_pred_xr        = 3;
base_dist.mean_vec = zeros(num_pred_xr, 1);
base_dist.scale    = 1 * eye(num_pred_xr);
base_dist.lambda   = 0.001;
base_dist.df       = num_pred_xr + 1;

prior_eta.v     = 2;
prior_eta.delta = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result   = cell(1, 1);
init     = cell(1, 10);
init_grp = 2;

time = (1:nknot)/(nknot+1);

for n_size=1:8
    rng(1132029);
    init{n_size}.weight      = ones(init_grp,1) / init_grp;
    init{n_size}.grp_cov_mat = zeros(num_pred_xr, num_pred_xr, init_grp);
    
    for i=1:init_grp
        init{n_size}.grp_cov_mat(:,:,i) = iwishrnd(base_dist.scale, base_dist.df);
    end
    
    init{n_size}.grp_mean = mvnrnd(base_dist.mean_vec, init{n_size}.grp_cov_mat /...
                           base_dist.lambda, init_grp)';
    init{n_size}.eta      = 1;
    init{n_size}.knot     = time;

    grp    = zeros(50+(n_size-1)*50,1);
    beta_i = zeros(num_pred_xr, 50+(n_size-1)*50);
    
    for i=1:(50+(n_size-1)*50)
        grp(i)      = randsample(init_grp, 1, true, init{n_size}.weight);
        beta_i(:,i) = mvnrnd(init{n_size}.grp_mean(:,grp(i)), ...
                      init{n_size}.grp_cov_mat(:,:,grp(i)));
    end
    
    init{n_size}.grp    = grp;
    init{n_size}.beta_i = beta_i;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preparing data and running MCMC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obs_time = df_random{nrep, nsize}(:,2);
ID       = df_random{nrep, nsize}(:,1);
y        = df_random{nrep, nsize}(:,3);
xf       = ones(length(ID),1);
xr       = df_random{nrep, nsize}(:,2);
        
for i=1:nknot
    tmp = df_random{nrep, nsize}(:,2) - time(i);
    negative = tmp < 0;
    tmp(negative) = 0;
    xr = [xr, tmp];
end
for i=1:(size(xr,2)-1)
    xr(:,i) = xr(:,i) - xr(:,(i+1));
end

same_indv = cell(ID(end),1);
for i=1:ID(end)
    same_indv{i} = find(ID == i);
end
Ti = cellfun('length', same_indv);
        
[result{1}.alpha, result{1}.sigma2, result{1}.beta_i, result{1}.r_int, ...
    result{1}.var_reff, result{1}.grp_mean, result{1}.grp_cov_mat, ...
    result{1}.knot, result{1}.grp, result{1}.active_grp, ...
    result{1}.weight, result{1}.eta] = mcmc(n_mcmc, ID, y, xf, xr, ...
        obs_time, same_indv, Ti, base_dist, prior_alpha, prior_var, ...
        prior_rint, prior_eta, init{nsize}, variable_knot);
    
est_grp_label = result{1}.grp(:,((burnin+thin):thin:n_mcmc));
