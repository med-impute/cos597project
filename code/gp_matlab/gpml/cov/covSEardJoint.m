function K = covSEardJoint(hyp, x, z, ii)
% Revised Squared Exponenetial Kernel for inter-patient modeling
% Squared Exponential covariance function with Automatic Relevance Detemination
% (ARD) distance measure. The covariance function is parameterized as:
% 
% Check kernel_rev.tex for reference
%
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
%% New input configuration of for joint modeling
% in the modified version, x = N x (V + 3), where
% N: total number of points (across individual, covariate, and time)
% 1: PAN (ID) of individual
% 2: Feature index
% 3: Timestamp (unit: accumulated hours)
% 4 to (V+3): V-dim demographic data
%
%% Get mode
if nargin<2, K = '1 + 2*D + 1/2*D*(D-1)'; return; end % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

%% Parse Data
disp('covSEardJoint: parsing data...')
[N, V] = size(x);
V = V - 3;
disp(['covSEardJoint: current # of dimension for inter-individual kernel: V = ', num2str(V)])

train_pid = x(:, 1);
train_pnum = length(unique(train_pid))

train_fid = x(:, 2);
train_fnum = 3 %length(unique(train_fid))

train_time = x(:, 3);
train_demogr = x(:, 4:(V+3))

if(~dg && ~xeqz)                                                     
    % get test cases
    [Ns, Vs] = size(z);
    Vs = Vs - 3;

    if(Vs ~= V)
        error('covSEardJoint: inconsistent number of dimenstion between training and testing cases')
    end
    test_pid = z(:, 1);
    test_pnum = length(unique(test_pid));

    test_fid = z(:, 2);
    test_fnum = 3 % length(unique(test_fid));

    test_time = z(:, 3);
    test_demogr = z(:, 4:(V+3));
end

%% Parse hyperparameter
disp('covSEardJoint: parsing hyperparam...')
lambda_ii = exp(hyp(1))
theta = exp(2*hyp(2:(train_fnum+1)))
lambda_j = exp(hyp((train_fnum+2):(2*train_fnum+1)))
lambda_jj = exp(squareform(hyp((2*train_fnum+2):end)))

%%
disp('covSEardJoint: computing kernels...')
% precompute squared distances 
if dg          
    disp('covSEardJoint: computing self cov. vector of testing cases...')
    % vector kxx/self-cov. vector of testing cases
    % (not sure)
    K = zeros(size(x, 1), 1);
else
  if xeqz        
    disp('covSEardJoint: computing self cov. of training cases...')
    % symmetric matrix Kxx/cov. matrix of training cases
    % initialization
%     K = zeros(N, N);
    Kc = ones(N, N);
    Kj = zeros(N, N);
    Kjj = zeros(N, N);
    Kii = zeros(N, N);
%     disp('start demogr dist...')
    if(V > 0)
        demogr_dist = squareform(pdist(train_demogr, 'correlation')); % optimize memory space later
    end
    
%     disp('covSEardJoint: start cross cov. term...')
    for i = 1:N
        for j = 1:N
            
            Kc(i, j) = 1/2; % theta(train_fid(j))/theta(train_fid(i));
            
            % if i == j: diag. case == 0
            if(i ~= j)
                if(train_pid(j) == train_pid(i))
                    % intra-individual case
                    if(train_fid(j) == train_fid(i))
                        % intra-covariate case
                        if(train_time(j) == train_time(i))
                            error(['dulplicate time stamps for patient: ', num2str(train_pid(j))])
                        else
                            Kj(i, j) = ((train_time(j) - train_time(i))/lambda_j(train_fid(j)))^2;
                        end
                    else
                        % inter-covariate case
                        if(train_time(j) == train_time(i))
                            % need to find out correlation for
                            % inter-covariate cases
                            Kjj(i, j) = ((0.4)/lambda_jj(train_fid(i), train_fid(j)))^2;                  
                        else
                            Kj(i, j) = (train_time(j) - train_time(i))^2/lambda_j(train_fid(i))/lambda_j(train_fid(j));
                            Kjj(i, j) = ((0.4)/lambda_jj(train_fid(i), train_fid(j)))^2;                            
                        end
                    end
                else
                    Kj(i, j) = 10^20;
                    Kjj(i, j) = 10^20;
                    % inter-individual case
                    if(train_fid(j) == train_fid(i))
                        % intra-covariate case
%                         Kii(i, j) = 10; % (demogr_dist(i, j)/lambda_ii)^2;
                        Kii(i, j) = (0.2/lambda_ii)^2;
                    else
                        % inter-covariate case
                        Kii(i, j) = 10^2; % near-zero correlation
                    end
                end
            end
        end
    end
    
  else
    disp('covSEardJoint: computing cross cov. of training and testing cases...')
    % cross covariances Kxz/cov. matrix of training and testing cases
    % initialization
%     K = zeros(N, Ns);
    Kc = zeros(N, Ns);
    Kj = zeros(N, Ns);
    Kjj = zeros(N, Ns);
    Kii = zeros(N, Ns);
    % inter-individual case
    for i = 1:N
        for j = 1:Ns
            Kc(i, j) = 1/2; %theta(test_fid(j))/theta(train_fid(i));
            if(test_pid(j) == train_pid(i))
                % intra-individual case
                if(test_fid(j) == train_fid(i))
                    % intra-covariate case
                    if(test_time(j) == train_time(i))
                        Kj(i, j) = 0; % optimize dulplicate cases later
                    else
                        Kj(i, j) = ((test_time(j) - train_time(i))/lambda_j(test_fid(j)))^2;
                    end
                else
                    % inter-covariate case
                    if(test_time(j) == train_time(i))
                        % need to find out correlation for
                        % inter-covariate cases
                        Kjj(i, j) = ((0.4)/lambda_jj(train_fid(i), test_fid(j)))^2;                  
                    else
                        Kj(i, j) = (test_time(j) - train_time(i))^2/lambda_j(train_fid(i))/lambda_j(test_fid(j));
                        Kjj(i, j) = ((0.4)/lambda_jj(train_fid(i), test_fid(j)))^2;
                    end
                end
            else
                Kj(i, j) = 10;
                Kjj(i, j) = 10;
                % inter-individual case
                if(test_fid(j) == train_fid(i))
                    % intra-covariate case
                    % curr_demogr_dist = corr(test_demogr(j)', train_demogr(i)');
                    % Kii(i, j) = (curr_demogr_dist/lambda_ii)^2;
                    Kii(i, j) = (0.2/lambda_ii)^2;
%                     Kii(i, j) = 10;
                else
                    % inter-covariate case
                    Kii(i, j) = 10^2; % near-zero correlation
                end
            end
        end
    end
  end
end

%% Calculate partial derivatives (based on i-th hyperparameter)
if ~dg                                                         
    Kvj = exp(-Kj/2);
    Kvjj = exp(-Kjj/2);
    Kvii = exp(-Kii/2);
    K = Kc.*(Kvj.*Kvjj);%  + Kvii);

    if(xeqz)
        [~,p] = chol(K);
        if(p > 0)
            disp('covSEardJoint: K is positive definite')
        else
            disp('covSEardJoint: K is not positive definite')
        end
        figure(1)
        subplot(2, 2, 1), imagesc(K), title('K'), colorbar
        subplot(2, 2, 2), imagesc(Kvj), title('Kvj'), colorbar
        subplot(2, 2, 3), imagesc(Kvjj), title('Kvjj'), colorbar
        subplot(2, 2, 4), imagesc(Kvii), title('Kvii'), colorbar
    else
        figure(2)
        subplot(2, 2, 1), imagesc(K), title('K'), colorbar
        subplot(2, 2, 2), imagesc(Kvj), title('Kvj'), colorbar
        subplot(2, 2, 3), imagesc(Kvjj), title('Kvjj'), colorbar
        subplot(2, 2, 4), imagesc(Kvii), title('Kvii'), colorbar
    end
    
end

% covariance
if(nargin > 3)                                      % derivatives
    error('covSEardJoint: methods to compute derivatives are not complete')
    if(ii == 1)
        K = 2*K;
    elseif(ii <= (2 + train_fnum + 1/2*train_fnum*(train_fnum-1))) % length scale parameters
        if dg                                                     
            % vector kxx/self-cov. vector of testing cases
            % (not sure)
            K = K*0;
        else
            if(ii == 2) % ii term
                K = Kc*(Kvii.*(Kii/(exp(hyp(2))^2)));
                
            elseif(ii <= (2 + train_fnum)) % j term
                K = Kc*(Kvj.*(Kj/(exp(hyp(ii))^2)).*Kvjj);
                
            else % jj term
                K = Kc*(Kvj.*Kvjj.*(Kjj/(exp(hyp(ii))^2)));                
            end
        end           
    else
        error('Unknown hyperparameter')
    end
end
