function t = subsasgn(t,s,b)
%SUBSASGN Subscripted assignement for cp_tensor.
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also CP_TENSOR.


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

switch s(1).type    
    case '.'
        switch s(1).subs
            case 'lambda'
		if length(s) == 1
		    t = cp_tensor(b, t.u);
		else
		    t.lambda = subsasgn(t.lambda, s(2:end), b);
		end
            case {'u','U'}
		if length(s) == 1
		    t = cp_tensor(t.lambda, b);
		else
		    tmpu = subsasgn(t.u, s(2:end), b);
		    t = cp_tensor(t.lambda, tmpu);
		end
            otherwise
                error(['No such field: ', s(1).subs]);
        end      
    case '()'
	error('Cannot change individual entries in a CP tensor.')
    case '{}'   
	new_s(1).type = '.';
	new_s(1).subs = 'u';
	new_s(2:length(s)+1) = s;
	t = subsasgn(t, new_s, b);
    otherwise
        error('Invalid subsasgn.');
end


