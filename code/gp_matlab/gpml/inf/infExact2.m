function [post nlZ dnlZ] = infExact(hyp, mean, cov, lik, x, y)

% Exact inference for a GP with Gaussian likelihood. Compute a parametrization
% of the posterior, the negative log marginal likelihood and its derivatives
% w.r.t. the hyperparameters. See also "help infMethods".
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2014-06-10.
%                                      File automatically generated using noweb.
%
% See also INFMETHODS.M.

if iscell(lik), likstr = lik{1}; else likstr = lik; end
if ~ischar(likstr), likstr = func2str(likstr); end
if ~strcmp(likstr,'likGauss')               % NOTE: no explicit call to likGauss
  error('Exact inference only possible with Gaussian likelihood');
end
 
[n, D] = size(x);
disp('infExact: handling cov')
K = feval(cov{:}, hyp.cov, x);                      % evaluate covariance matrix
disp('infExact: handling mean')
m = feval(mean{:}, hyp.mean, x);                          % evaluate mean vector

sn2 = exp(2*hyp.lik);                               % noise variance of likGauss
if sn2<1e-6                        % very tiny sn2 can lead to numerical trouble
  disp('infExact: tiny sn2')
  L = chol(K+sn2*eye(n)); sl =   1;   % Cholesky factor of covariance with noise
  pL = -solve_chol(L,eye(n));                            % L = -inv(K+inv(sW^2))
else
  disp(['infExact: using sn2 = ', num2str(sn2)])
  disp(['infExact: n = ', num2str(n)])
  M = (K/sn2+eye(n));
  figure(2), imagesc(M), colorbar, title('M')
  
  disp(['infExact: L = chol(M)'])
  L = chol(K/sn2+eye(n)); sl = sn2;                       % Cholesky factor of B

  disp(['infExact: pL = ', num2str(L)])
  pL = L;                                           % L = chol(eye(n)+sW*sW'.*K)
end
disp('infExact: solve_chol')
alpha = solve_chol(L,y-m)/sl;

post.alpha = alpha;                            % return the posterior parameters
post.sW = ones(n,1)/sqrt(sn2);                  % sqrt of noise precision vector
post.L = pL;

if nargout>1                               % do we want the marginal likelihood?
  nlZ = (y-m)'*alpha/2 + sum(log(diag(L))) + n*log(2*pi*sl)/2;   % -log marg lik
  if nargout>2                                         % do we want derivatives?
    disp('infExact: allocate derivatives')
    dnlZ = hyp;                                 % allocate space for derivatives
    Q = solve_chol(L,eye(n))/sl - alpha*alpha';     % precompute for convenience
    for i = 1:numel(hyp.cov)
      dnlZ.cov(i) = sum(sum(Q.*feval(cov{:}, hyp.cov, x, [], i)))/2;
    end
    dnlZ.lik = sn2*trace(Q);
    for i = 1:numel(hyp.mean), 
      dnlZ.mean(i) = -feval(mean{:}, hyp.mean, x, i)'*alpha;
    end
  end
end
