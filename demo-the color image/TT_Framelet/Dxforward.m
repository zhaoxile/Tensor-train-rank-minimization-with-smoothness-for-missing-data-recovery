function ux=Dxforward(u,BoundaryCondition)
%用于求D_x^+, 默认循环边界
% ux=imfilter(u,[0,-1,1],BoundaryCondition);
%示例：ux=Dxforward(u,'circular');
if nargin<2
    BoundaryCondition='circular';
end
switch BoundaryCondition
    case 'circular'
        ux=[diff(u,1,2) u(:,1)-u(:,end)];
    case 'symmetric'
        ux=[diff(u,1,2) zeros(size(u,1),1)];
    case 'zero'
        ux=[diff(u,1,2) -u(:,end)];
end
end