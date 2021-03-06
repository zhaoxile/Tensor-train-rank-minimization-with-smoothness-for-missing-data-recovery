function c = ttv(a,v,n)
%TTV Tensor times vector for cp_tensor.
%
%   B = TTV(A,V,N) computes the product of CP tensor A with a (column)
%   vector V.  The integer N specifies the dimension in A along which
%   V is multiplied.  If size(V) = [I,1], then A must have size(A,N) =
%   I.  Note that ndims(B) = ndims(A) - 1 because the N-th dimension
%   is removed.
%
%   B = TTV(A,U) computes the product of tensor A with a sequence of
%   vectors in the cell array U.  The products are computed
%   sequentially along all dimensions (or modes) of A. The cell array
%   U contains ndims(A) vectors.
%
%   B = TTV(A,U,DIMS) computes the sequence tensor-vector products
%   along the dimensions specified by DIMS.
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also TENSOR/TTV, CP_TENSOR.

%Brett W. Bader and Tamara G. Kolda, Released under SAND2004-5189,
%Sandia National Laboratories, 2004.  Please address questions or
%comments to: tgkolda@sandia.gov.  Terms of use: You are free to copy,
%distribute, display, and use this work, under the following
%conditions. (1) You must give the original authors credit. (2) You may
%not use or redistribute this work for commercial purposes. (3) You may
%not alter, transform, or build upon this work. (4) For any reuse or
%distribution, you must make clear to others the license terms of this
%work. (5) Any of these conditions can be waived if you get permission
%from the authors.


%%%%%%%%%%%%%%%%%%%%%%
%%% ERROR CHECKING %%%
%%%%%%%%%%%%%%%%%%%%%%

% Check the number of arguments
if (nargin < 2)
    error('TTV requires at least two arguments.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CHECK FOR CELL ARRAY %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isa(v,'cell')
    if (nargin < 3)
	n = 1:ndims(a);
    end
    c = ttvs(a,v,n);
    return;
end

%%%%%%%%%%%%%%%%%%%%%%
%%% ERROR CHECKING %%%
%%%%%%%%%%%%%%%%%%%%%%

% Check the second argument
if (ndims(v) ~= 2) || (size(v,2) ~= 1) || (size(v,1) == 1) 
    error('V must be a column vector.');
end

% 3rd argument
if (nargin < 3) || ~isequal(size(n),[1 1]) || (n < 0) || (n > ndims(a))
    error('Dimension N must be between 1 and NDIMS(A).');
end

% Check that sizes match!
if size(a,n) ~= size(v,1)
    error('Nth dimension of A is not equal to length of V.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPUTE THE PRODUCT %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute inner product with the nth mode for each of the factors 
tmp = a.u{n}'*v;
newLambda = a.lambda .* tmp;

% Check if result is a scalar
if ndims(a) == 1
    c = sum(newLambda);
    return;
end

% Otherwise, store as new CP tensor
c = cp_tensor(newLambda, {a.u{1:n-1}, a.u{n+1:end}});

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c = ttvs(a,v,dims)

% Check validity of DIMS and then check for "minus" case
if (max(abs(dims)) > ndims(a))
    error('An entry in DIMS exceeds order of A.');
elseif (max(dims) < 0)
    dims = setdiff(1:ndims(a), -dims);
end
    
% Check validity of parameters passed to TTVS
n = length(dims);
if (n > ndims(a)) || (n > length(v))
    error('DIMS is too long.');
elseif (n < length(v)) && (length(v) < ndims(a))
    error('If length(DIMS) < length(V), then length(V) must equal ndims(A).');
elseif (length(v) > ndims(a))
    error('Length of V greater than order of A.');
end

% Reorder dims from largest to smallest
[sdims,sidx] = sort(-dims);
sdims = -sdims;

% Check sizes of U and DIMS to determine version to use.
if (n == length(v))
    vidx = sidx;   % index V by sorted order 
else 
    vidx = sdims;  % index V by (sorted) dimension 
end

% Calculate individual products
c = ttv(a, v{vidx(1)}, sdims(1));
for i = 2 : n
    c = ttv(c, v{vidx(i)}, sdims(i));
end
