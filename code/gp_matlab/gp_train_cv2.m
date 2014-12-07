%% Setup hold-out set and "corresponding training set"
mse_score = NaN*zeros(test_pnum, test_fnum);
corr_score = NaN*zeros(test_pnum, test_fnum);

% Leave-One-Out Cross Validation
for ti = 1:train_pnum
    curr_inter = (train_x(:, 1) ~= train_id(ti));

    for tj = 1:fnum
        flag = 1;

        m2 = [];
        s2 = [];
        x2 = train_x((train_x(:, 1) == train_id(ti)) & (train_x(:, 2) == tj), 3);
        y2 = train_y((train_x(:, 1) == train_id(ti)) & (train_x(:, 2) == tj));
        curr_resample =  train_z(...
                            (train_z(:, 1) == train_id(ti))...
                          & (train_z(:, 2) == tj), :);

        for tt = 1:length(curr_resample)
            % select reasonable training set
            % add observed sample from the same patient
            curr_intra =    (train_x(:, 1) == train_id(ti))...
                          & (train_x(:, 3) < curr_resample(tt, 3));

            curr_x = [train_x(curr_intra, :); train_x(curr_inter, :)];
            curr_y = [train_y(curr_intra, :); train_y(curr_inter, :)];

            curr_z = curr_resample(tt, :);

            % GP training
            try
                hyp = minimize(hyp, @gp, -iter_num, @infExact, meanfunc, covfunc, likfunc, curr_x, curr_y);
                [m s] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, curr_x, curr_y, curr_z);

                % try FITC
                % [m s] = gp(hyp, @infFITC, meanfunc, covfuncF, likfunc, x, y, z);
                % display the results
                % record mse/corr
                m2 = [m2; m];
                s2 = [s2; s];

            catch
                disp('error in this round...record and skipping...');
                failure_set = [failure_set; cov_vec'];
                flag = 0;
                break;
            end
        end

        % successfully iterate all time points,
        % calculate MSE/correlation 
        if(flag == 1)
            curr_check_id = train_checkpt{ti, tj};
            if(~isempty(y2))
                if(length(y2) ~= length(curr_check_id))
                    error('mismatch in testing cases')
                end
                curr_mse = sum((y2 - m2(curr_check_id)).^2)/length(curr_check_id);
                curr_mse = sqrt(curr_mse);
%                                         curr_mse = curr_mse*(global_std(tj)^2);            
                if(length(curr_check_id) > 1)
                    curr_corr_score = corrcoef(y2, m2(curr_check_id));
                    curr_corr_score = curr_corr_score(1, 2);
                else
                    curr_corr_score = NaN;
                end
                mse_score(ti, tj) = curr_mse;
                corr_score(ti, tj) = curr_corr;
            end

            % output regression figure for reference
            if(write_fig == 1)
                z2 = curr_resample(:, 3);

                f = [((m2+2*sqrt(s2))*post_scale(tj) + global_mean(tj));...
                    flipdim(((m2-2*sqrt(s2))*post_scale(tj) + global_mean(tj)),1)];
                a = fill([z2; flipdim(z2,1)], f, [7 7 7]/8);
                hold on;
                set(a,'EdgeColor','none');

                plot(z2, (m2.*post_scale(tj) + post_offset(tj)), 'LineWidth', 2); 
                plot(x2, (y2.*post_scale(tj) + post_offset(tj)), '+', 'MarkerSize', 12);
                grid on

                xlabel('input, elspased hour')
                ylabel(['output, ', new_cat{select_fidx(tj)}])
                title(['Joint GP Regression Result: patient ', num2str(train_id(ti)),...
                    ', var ', new_cat{select_fidx(tj)}])
                hold off;
                print(['figure_joint/', 'gp_joint_PAN_', num2str(train_id(ti)),...
                    '_VAR_', num2str(tj), '_', new_cat{select_fidx(tj)}], '-djpeg'); 
            end

        end

    end
end

%% record two best_set through averaging MSE and corr
% evaluate through MSE 
% (average of each feature and product them together)
curr_avg_mse = mean(mse_score(~isnan(mse_score)));
if(curr_avg_mse < curr_best_mse)
    curr_best_mse = curr_avg_mse;
    best_set_mse = exp(hyp.cov);
end
% evaluate through correlation
% (average of each feature and product them together)
curr_avg_corr = mean(corr_score(~isnan(corr_score)));
if(curr_avg_corr > curr_best_corr)
    curr_best_corr = curr_avg_corr;
    best_set_corr = exp(hyp.cov);
end
 