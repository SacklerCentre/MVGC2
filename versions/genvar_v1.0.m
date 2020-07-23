%% genvar
%
% Generate VAR time series data
%
% <matlab:open('genvar.m') code>
%
%% Syntax
%
%     [X,E] = genvar(A,E,trunc)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     A          VAR coefficients matrix
%     E          residuals (single-trial) time series
%     trunc      number of initial observations to truncate (default: 0)
%
% _output_
%
%     X          VAR (single-trial) time series with residuals E
%     E          possibly truncated residuals
%
%% Description
%
% Generates a vector autoregression (VAR) time series with coefficients |A| and
% residuals |E|; implements:
%
% <<eq_varnn.png>>
%
% Initial values for the outputs |X| are set to the residuals |E|. Optionally
% truncate the first |trunc| observations.
%
% *_Note:_* If available will use the C routine <matlab:open('C/genvar_mex.c')
% |genvar_mex.c|> to perform the actual computation; otherwise a (slower)
% Matlab-coded routine is used. The existence of a |genvar| mex file for the
% current platform is checked for in the <startup.html |startup|> script.
%
%% See also
%
% <var_to_tsdata.html |var_to_tsdata|> |
% <startup.html |startup|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [X,E] = genvar_v1.0(A,E,trunc)

global have_genvar_mex;

assert(isnumeric(E) && isreal(E) && ismatrix(E), 'Residuals must be a real matrix')

[n,m] = size(E);

if nargin < 3 || isempty(trunc)
    trunc = 0;
else
    assert(trunc >= 0 && trunc < m,'bad truncation');
end

if isempty(A) % model order p = 0, just return residuals (v2.0)

    X = E;

else

    assert(isnumeric(A) && isreal(A) && (ismatrix(A) || ndims(A) == 3), 'VAR coefficients must be a real matrix or 3D array')

    [n1,n2,p] = size(A);

    assert(n1 == n2, 'VAR coefficients blocks not square');
    assert(n1 == n,  'Residuals covariance matrix doesn''t match VAR coefficients matrix');

    if have_genvar_mex
        X = genvar_mex(A,E);
    else
        X = E;
        for t = p+1:m
            for k = 1:p
                X(:,t) = X(:,t) + A(:,:,k)*X(:,t-k);
            end
        end
    end

end

if trunc > 0
    X = X(:,trunc+1:m);
    if nargout > 1
        E = E(:,trunc+1:m);
    end
end
