%% ME 6000 Numerical Analysis - Assignment for Numerical Errors
%  Taylor series of exp(x):  e^x = sum_{k=0}^{N} x^k / k!
%
%  For x = -10 and x = +10, N = 10, 30, 50, this script performs:
%    (1) double precision, increasing k   (k = 0,1,2,...,N)
%    (2) double precision, decreasing k   (k = N,N-1,...,0)
%    (3) single precision, increasing k
%    (4) single precision, decreasing k
%  and compares every result with the exact value exp(x).
%
%  NOTE on how the terms are formed
%  --------------------------------
%  Computing x^k/factorial(k) literally overflows in single precision:
%  for x = -10, N = 50 we would need 10^50 and 50! ~ 3.04e64, while the
%  largest single-precision number is only ~3.40e38  ->  Inf/Inf = NaN.
%  We therefore build the terms with the stable recurrence
%        t_0 = 1,      t_k = t_{k-1} * x / k
%  which is mathematically identical but keeps every intermediate value
%  inside the representable range (max |t_k| = 2755.7 for x = -10).
%  The rounding behaviour we want to study is unaffected.

clear; clc; format long e

xList = [-10, 10];
NList = [10, 30, 50];

for x = xList
    exact = exp(x);                       % reference ("exact") value
    fprintf('\n==================================================================\n');
    fprintf(' x = %+d      exact e^x = %.15e\n', x, exact);
    fprintf('==================================================================\n');
    fprintf('%5s  %-10s %-12s %24s %14s %14s\n', ...
            'N','precision','order','summation','abs. error','rel. error [%]');

    for N = NList
        % ---- the four required cases -------------------------------------
        Sdi = taylorSum(x, N, 'double', 'increasing');   % (1)
        Sdd = taylorSum(x, N, 'double', 'decreasing');   % (2)
        Ssi = taylorSum(x, N, 'single', 'increasing');   % (3)
        Ssd = taylorSum(x, N, 'single', 'decreasing');   % (4)

        printRow(N,'double','increasing',Sdi,exact);
        printRow(N,'double','decreasing',Sdd,exact);
        printRow(N,'single','increasing',Ssi,exact);
        printRow(N,'single','decreasing',Ssd,exact);
        fprintf('%s\n', repmat('-',1,90));
    end
end

%% ---------------------------------------------------------------------
function S = taylorSum(x, N, prec, order)
% TAYLORSUM  Partial Taylor sum of exp(x) up to k = N.
%   prec  : 'double' or 'single'   -> precision of EVERY operation
%   order : 'increasing' or 'decreasing' order of accumulation

    cast = str2func(prec);            % @double or @single
    xc   = cast(x);

    % --- build the terms t_k = x^k/k! with the stable recurrence --------
    t    = zeros(1, N+1, prec);
    t(1) = cast(1);                                   % k = 0
    for k = 1:N
        t(k+1) = t(k) * xc / cast(k);                 % t_k = t_{k-1}*x/k
    end

    % --- accumulate in the requested order ------------------------------
    S = cast(0);
    if strcmp(order,'increasing')
        idx = 1:N+1;                  % k = 0,1,...,N   (smallest first)
    else
        idx = N+1:-1:1;               % k = N,...,1,0   (smallest last)
    end
    for i = idx
        S = S + t(i);                 % stays in the chosen precision
    end
end

%% ---------------------------------------------------------------------
function printRow(N, prec, order, S, exact)
    absErr = abs(double(S) - exact);
    relErr = absErr/abs(exact)*100;
    fprintf('%5d  %-10s %-12s %24.15g %14.4e %14.6g\n', ...
            N, prec, order, double(S), absErr, relErr);
end
