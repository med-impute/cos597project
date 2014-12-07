%% Testing flow: 
% (1) pick 100 sepsis patients as training set
% (2) perform cross validation for hyper parameter picking
% (3) pick another 100 sepsis patients as testing set
% (4) perform "on-line" testing and calculate the score (both MSE and correlation)
load('/memex/lifangc/COS597D_Project/share/raw_data_5000_patient.mat')

LIB = './';
GPML = '/gpml';
path(path, [LIB, GPML])
path(path, [LIB, GPML, '/cov'])
path(path, [LIB, GPML, '/doc'])
path(path, [LIB, GPML, '/inf'])
path(path, [LIB, GPML, '/lik'])
path(path, [LIB, GPML, '/mean'])
path(path, [LIB, GPML, '/util'])

%% Parameter Setting
write_fig = 1;

train_pnum = 5;
test_pnum = 5;

select_fidx = [1, 2, 3]; % access variable 'new_cat' for detailed feature info.
fnum = length(select_fidx);
fmat = FINAL_FEATURE_MAT(:, select_fidx);

sepsis_idx = find(FINAL_SEPSIS_FLAG == 1);
s = rng;
select_idx = randperm(length(sepsis_idx), length(sepsis_idx));
rng(s);

%% De-mean for the feature matrix first (use "training set")
% i: individual, j: covariate, t: time
for j = 1:length(select_fidx)
    global_mean(j) = ...
        mean(fmat(~isnan(fmat(:, j)), j));
    global_std(j) = ...
        std(fmat(~isnan(fmat(:, j)), j));

    fmat(:, j) = fmat(:, j) - global_mean(j);
%     fmat(:, j) = fmat(:, j)/global_std(j);
    
%     mean(fmat(~isnan(fmat(:, j)), j))
%     std(fmat(~isnan(fmat(:, j)), j))
    
end

%% Run scripts
gp_load_data
gp_train_joint2
% gp_test_joint       
         


