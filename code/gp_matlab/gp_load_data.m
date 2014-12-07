%% Pick training and testing data
% x: (patient id, feature id, time) + demographic data (TODO)
% y: regression value (true observed value)
% z: regression time value (smoothing time)

%% Get training data
train_x = [];
train_y = [];
train_z = [];
train_id = [];
train_checkpt = cell(train_pnum, fnum);

count_p = 1;
count_id = 1;
while(count_p <= train_pnum)
    curr_id = VISIT_ID(select_idx(count_id));
    block_id = find(FINAL_ID == curr_id);
    
    if(~isempty(block_id))
        train_id = [train_id, curr_id];
        tmp_hour = (FINAL_TIME(block_id) - FINAL_TIME(block_id(1)))/3600;
        
        % get regression time
        resample_t = [];
        for t = 1:length(block_id)
            if(t == length(block_id))
                resample_t = [resample_t; tmp_hour(t)];
                for j = 1:fnum
                    if(~isnan(fmat(block_id(t), j)))
                        train_checkpt{count_p, j} = [train_checkpt{count_p, j}, length(resample_t)];
                    end
                end
            else
                interp_time = linspace(tmp_hour(t),...
                                       tmp_hour(t+1),...
                                       floor(tmp_hour(t+1) - tmp_hour(t)));
                if(isempty(interp_time))
                    interp_time = [tmp_hour(t), tmp_hour(t+1)];
                end

                for j = 1:fnum
                    if(~isnan(fmat(block_id(t), j)))
                        train_checkpt{count_p, j} = [train_checkpt{count_p, j}, length(resample_t) + 1];
                    end
                end
                resample_t = [resample_t; interp_time(1:(end-1))'];
            end
        end
        
        % get data
        for j = 1:fnum
            valid_id = find(~isnan(fmat(block_id, j)));
            train_x = [...
                        train_x;...
                        curr_id*ones(length(valid_id), 1),...
                        j*ones(length(valid_id), 1),...
                        tmp_hour(valid_id)...
                        ]; % trainining cases
            train_y = [...
                        train_y;...
                        fmat(block_id(valid_id), j)...
                        ]; % training values
            train_z = [...
                        train_z;...
                        curr_id*ones(length(resample_t), 1),...
                        j*ones(length(resample_t), 1),...
                        resample_t...
                        ]; % testing cases
        end
          
        % update counter
        count_p = count_p + 1;
    end
    count_id = count_id + 1;
end

%% Get testing data
test_x = [];
test_y = [];
test_z = [];
test_id = [];
test_checkpt = cell(test_pnum, fnum);

count_p = 1;
while(count_p <= test_pnum)
    curr_id = VISIT_ID(select_idx(count_id));
    block_id = find(FINAL_ID == curr_id);
    
    if(~isempty(block_id))
        test_id = [test_id, curr_id];
        tmp_hour = (FINAL_TIME(block_id) - FINAL_TIME(block_id(1)))/3600;
        
        % get regression time
        resample_t = [];
        for t = 1:length(block_id)
            if(t == length(block_id))
                resample_t = [resample_t; tmp_hour(t)];
                for j = 1:fnum
                    if(~isnan(fmat(block_id(t), j)))
                        test_checkpt{count_p, j} = [test_checkpt{count_p, j}, length(resample_t)];
                    end
                end
            else
                interp_time = linspace(tmp_hour(t),...
                                       tmp_hour(t+1),...
                                       floor(tmp_hour(t+1) - tmp_hour(t)));
                if(isempty(interp_time))
                    interp_time = [tmp_hour(t), tmp_hour(t+1)];
                end

                for j = 1:fnum
                    if(~isnan(fmat(block_id(t), j)))
                        test_checkpt{count_p, j} = [test_checkpt{count_p, j}, length(resample_t) + 1];
                    end
                end
                resample_t = [resample_t; interp_time(1:(end-1))'];
            end
        end
        
        % get data
        for j = 1:fnum
            valid_id = find(~isnan(fmat(block_id, j)));
            test_x = [...
                        test_x;...
                        curr_id*ones(length(valid_id), 1),...
                        j*ones(length(valid_id), 1),...
                        tmp_hour(valid_id)...
                        ]; % trainining cases
            test_y = [...
                        test_y;...
                        fmat(block_id(valid_id), j)...
                        ]; % training values
            test_z = [...
                        train_z;...
                        curr_id*ones(length(resample_t), 1),...
                        j*ones(length(resample_t), 1),...
                        resample_t...
                        ]; % testing cases
        end
          
        % update counter
        count_p = count_p + 1;
    end
    count_id = count_id + 1;
end

