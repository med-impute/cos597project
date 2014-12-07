%% Joint GP training and parameter selection using cross validation
% hyp = [ 
%         log(lambda_ii)  // 1 : class-wise individual length scale (1 unique coef. so far)
%         log(theta_1)    // 1 + (1 to D): signal variance for each covariate
%         log(theta_2)  
%          ...
%         log(theta_D)    
%         log(lambda_1)   // 1 + D + (1 to D): temporal length scale for each covariate
%         log(lambda_2)  
%          ...
%         log(lambda_D)    
%         log(lambda_12)   // 1 + 2*D + (1 to 1/2*D*(D-1)):
%         log(lambda_13)   // pairwise length scale for each pair of covariate
%          ...
%         log(lambda_D-1D) // Total # of hyperparameter now: 1 + 2*D + 1/2*D*(D-1)
%                         ]
clear hyp
lik_level = [0.1, 1, 10, 100, 1000]; % noise level need to adjust later as well
j_level = [0.1, 1, 10, 100, 1000]; % temporal length scale to test
jj_level = [0.1, 1, 10, 100, 1000]; % 0.1:0.2:0.9; % inter-covariate level to test
ii_level = [0.1, 1, 10, 100, 1000]; % 0.1:0.2:0.9; % inter-individual level to test

post_scale = ones(length(select_fidx), 1);
post_offset = zeros(length(select_fidx), 1);
% post_scale = global_std;
post_offset = global_mean;

%% Fixed function and value
iter_num = 1000;
meanfunc = @meanZero;
likfunc = @likGauss; 
 
covfunc = @covSEardJoint;
% try infFITC
% n = size(train_x, 1);
% nu = fix(n/2); iu = randperm(n); iu = iu(1:nu); 
% u = train_x(iu,:);
% covfuncF = {@covFITC, {covfunc}, u};

%% Setup hyp.cov
curr_best_mse = Inf;
curr_best_corr = 0;
best_set_mse = [];
best_set_corr = [];
failure_set = [];
for ii = ii_level
    for j1 = j_level
        for j2 = j_level
            for j3 = j_level
                for jj = jj_level
                    for lik = lik_level
                        tic
                        %% Setup hyperparamenter for covariance function
                        hyp.lik = log(lik);
                        cov_vec = [...
                                    ii;...
                                    (global_std)';...
                                    j1; j2; j3;...
                                    jj*ones(1/2*fnum*(fnum-1), 1);...
                                    ];
                        hyp.cov = log(cov_vec);
                        
                        try
                            hyp = minimize(hyp, @gp, -iter_num, @infExact, meanfunc, covfunc, likfunc, train_x, train_y);
                            [m s] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, train_x, train_y, train_x);
                            
                                                    
                            gp_train_cv
                            
                        catch
                            lasterr
                            disp('error in this round...record and skipping...');
                            failure_set = [failure_set; cov_vec'];
                            flag = 0;
                            break;                            
                        end
                        
                        toc
                    end
                end
            end
        end
    end
end

% if(write_fig == 1)
% mse_score = zeros(test_pnum, test_fnum);
% corr_score = zeros(test_pnum, test_fnum);
% for testp = 1:test_pnum
% for testf = 1:test_fnum
% figure(1000 + testf)
% m2 = m((z(:, 2) == testf) & (z(:, 1) == VISIT_ID(testp)));
% s2 = s((z(:, 2) == testf) & (z(:, 1) == VISIT_ID(testp)));
% z2 = z((z(:, 2) == testf) & (z(:, 1) == VISIT_ID(testp)), 3);
% x2 = x((x(:, 2) == testf) & (x(:, 1) == VISIT_ID(testp)), 3);
% y2 = y((x(:, 2) == testf) & (x(:, 1) == VISIT_ID(testp)));
% 
% curr_check_id = check_table{testp, testf};
% if(~isempty(y2))
%     if(length(y2) ~= length(curr_check_id))
%         error('mismatch in testing cases')
%     end
%     curr_mse = sum((y2 - m2(curr_check_id)).^2)/length(curr_check_id);
%     curr_mse = curr_mse*(global_std(testf)^2);            
%     if(length(curr_check_id) > 1)
%         curr_corr_score = corrcoef(y2, m2(curr_check_id));
%         curr_corr_score = curr_corr_score(1, 2);
%     else
%         curr_corr_score = 0;
%     end
% else
%     curr_mse = Inf;
%     curr_corr_score = 0;
% end
% mse_score(testp, testf) = curr_mse;
% corr_score(testp, testf) = curr_corr_score;
% 
% f = [((m2+2*sqrt(s2))*post_scale(testf) + global_mean(testf));...
%     flipdim(((m2-2*sqrt(s2))*post_scale(testf) + global_mean(testf)),1)];
% a = fill([z2; flipdim(z2,1)], f, [7 7 7]/8); hold on;
% set(a,'EdgeColor','none');
% %         set(a,'facealpha',0.6);
% 
% plot(z2, (m2.*post_scale(testf) + post_offset(testf)), 'LineWidth', 2); 
% plot(x2, (y2.*post_scale(testf) + post_offset(testf)), '+', 'MarkerSize', 12);
% grid on
% %     xlim([min(z), max(x)])
% xlabel('input, elspased hour')
% ylabel(['output, ', new_cat{testf}])
% title(['Joint GP Regression Result: MSE = ', num2str(curr_mse),...
%         ', Correlation Score = ', num2str(curr_corr_score)])
% hold off;
% print([DIR_FIGURE, 'gp_joint_PAN_', num2str(VISIT_ID(testp)),...
%     '_VAR_', num2str(testf), '_', new_cat{testf}], '-depsc'); 
% end
% end
% end


%% Plot optimal training results
hyp.cov = log([0.1;...
               (global_std)';...
               4; 4; 4;...
               0.5*ones(1/2*fnum*(fnum-1), 1);...
            ]);
        
        
        
