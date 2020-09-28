function [R] = quasimvnrnd(mu,sigma,n,seq)
%Calculates Quasi-random sample of a multivariate normal distribution with
%expectation vector mu and covariance matrix sigma based on the
%quasi-random sequence specified by "seq"
%Input: mu - (dx1) vector of expectation value
%       sigma -(dxd) covariance matrix
%       n - size of quasi-random sample
%       seq - choose "halton" to use Halton sequence, "sobol" for Sobol
%             sequence and "latin" for Latin hypercube sequence
%  Output:  R - (nxd)matrix of n realisations

d=size(mu,2);
if strcmp(seq,'halton')
p = haltonset(d,'Skip',1e3,'Leap',1e2);
p = scramble(p,'RR2');
X = net(p,n);
elseif strcmp(seq,'sobol')
p = sobolset(d,'Skip',1e3,'Leap',1e2);
p = scramble(p,'MatousekAffineOwen');
X = net(p,n);
elseif strcmp(seq,'latin')
X= lhsdesign(n,d);
else
    disp('Please enter valid quasi-random sequence, choose "halton" to use Halton sequence, "sobol" for Sobol sequence and "latin" for Latin hypercube sequence.');
end
    
 
for i=1:d
    R(:,i)=norminv(X(:,i),mu(i),sqrt(sigma(i,i)));
end

end

