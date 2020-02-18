function uy=Dybackward(u,BoundaryCondition)
%用于求D_y^-, 默认循环边界
% uy=imfilter(u,[-1;1;0],BoundaryCondition);
%示例：uy=Dybackward(u,'circular');
if nargin<2
    BoundaryCondition='circular';
end
switch BoundaryCondition
    case 'circular'
        uy=[u(1,:)-u(end,:);diff(u,1,1)];
    case 'symmetric'
        uy=[zeros(1,size(u,2));diff(u,1,1)];
    case 'zero'
        uy=[u(1,:);diff(u,1,1)];
end
end