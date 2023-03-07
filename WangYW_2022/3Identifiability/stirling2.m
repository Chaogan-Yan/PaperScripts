function S = stirling2(m, n, c)
%STIRLING2 returns Stirling number of the 2nd kind
%
%BACKGROUND Stirling numbers of the 2nd kind show how many ways in which a
%set of m elements can be clustered to form n groups. This algorithm uses a
%recursive formula that is exact up to intmax('uint64'). Any Stirling
%number under this will be represented exactly and any number larger than
%this but less than realmax will be approximated using double precision
%(16-digits). This will solve Stirling numbers exactly for S(26,n), for any
%n and S(>26,n) for select n. Solutions are available for S(<220,n) for all
%n, and S(>=220,n) for select n. An error will indicate overflow. If the S
%is uint64 it is exact. If S is double, it may be an approximation good to
%16 digits.
%
%USAGE
%   S = Stirling( m ) 
%   returns S, a vector containing the number of ways in which a set of m
%   elements can be clustered to form 0,1,2,...,m groups
%
%   S = Stirling( m, n )
%   returns scalar S, the number of ways in which a set of m
%   elements can be clustered to form m groups
%
%   S = Stirling( m, n, 'double'). 
%   This convention is used internally, if it is empirically determined
%   that S > intmax('unit64'). You can use it directly to force double
%   precision.


if nargin < 2 || isempty(n)
    n = 0:m;
end
if nargin < 3 
    c = 'uint64';
end
S = diag(ones(m+1,1,c) );       % table of solutions S(m,n)
    
for r = 2:m+1
    for c = 2:r-1
        S(r,c) = S(r-1, c-1) + (c-1)*S(r-1,c);
    end
end

if nargin > 1
    S = S(end,n+1);
else
    S = S(end,:);
end

if isa(S,'double') && any(S > realmax('double'))
    error('overflow');
elseif isa(S,'uint64') && any(S >= intmax('uint64'))
    S = stirling2(m,n,'double');
end

S = S(:);

