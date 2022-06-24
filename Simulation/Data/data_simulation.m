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

% Script for generating simulated data
df_fixed  = cell(100, 10);
df_random = cell(100, 10);

for nsize=1:10
    for iter=1:100

        rng(1132000+iter)
        P        = 50+(nsize-1)*50;
        weight   = [0.25; 0.25; 0.25; 0.25];
        sd_error = sqrt(0.15);
        t_min    = 10;
        t_max    = 20;
        index    = randsample(1:4, P, true, weight);
        P_grp    = tabulate(index);
        P_grp    = P_grp(:,2);
        Ti       = randsample(t_min:t_max, P, true)';
        ID       = repelem(1:P, 1, Ti)';

        reff_mean = 0.75;
        reff_var  = 0.5;
        reff      = randn(ID(end), 1) * sqrt(reff_var) + reff_mean;
        error     = randn(length(ID),1) * sd_error;

        grp_mean    = [-3.0, -7.5, -3.0,  4.0;
                       -3.0, -5.0, -1.0,  1.0;
                       -3.0,  0.0,  3.0, -3.0];
        grp_cov_mat = repmat(0.2 * corrcov(iwishrnd([1 -0.5 0; -0.5 1 -0.5; 0 -0.5, 1], 100)), [1 1 length(weight)]);

        beta_i = zeros(size(grp_mean, 1), P);
        beta_i(:,1:P_grp(1))                            = mvnrnd(grp_mean(:,1), grp_cov_mat(:,:,1),P_grp(1))';
        beta_i(:,(P_grp(1)+1):sum(P_grp(1:2)))          = mvnrnd(grp_mean(:,2), grp_cov_mat(:,:,2),P_grp(2))';
        beta_i(:,(sum(P_grp(1:2))+1):(sum(P_grp(1:3)))) = mvnrnd(grp_mean(:,3), grp_cov_mat(:,:,3),P_grp(3))';
        beta_i(:,(sum(P_grp(1:3))+1):end)               = mvnrnd(grp_mean(:,4), grp_cov_mat(:,:,4),P_grp(4))';

        same_indv = cell(ID(end),1);
        for i=1:ID(end)
            same_indv{i} = find(ID == i);
        end

        for knot_val=1:2
            xr = [];
            for i=1:P
                xr = [xr; sort(randsample(1:365, Ti(i), false))'];
            end
            xr          = xr / 365;
            obs_time    = xr;
            random_knot = knot_val;
            knot        = zeros(P, 2);
            if random_knot == 1
                for i=1:P
                    knot(i,:) = [rand(1) * 0.5, rand(1) * 0.5 + 0.5];
                end
            else
                knot(:,1) = 1/3;
                knot(:,2) = 2/3;
            end
            nknot = size(knot, 2);
            for i=1:nknot
                tmp           = obs_time - repelem(knot(:,i), Ti, 1);
                negative      = tmp < 0;
                tmp(negative) = 0;
                xr            = [xr, tmp];
            end
            for i=1:(size(xr,2)-1)
                xr(:,i) = xr(:,i) - xr(:,(i+1));
            end
            y = [];
            for i=1:P
                y = [y; reff(i)+xr(same_indv{i},:) * beta_i(:,i) + error(same_indv{i})];
            end
            data = [ID, obs_time, y];
            if random_knot == 1
                df_random{iter, nsize} = data;
            else
                df_fixed{iter, nsize} = data;
            end
        end
    end
end

save("df_fixed.mat", "df_fixed");
save("df_random.mat", "df_random");
