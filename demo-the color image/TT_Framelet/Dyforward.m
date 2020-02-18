function uy=Dyforward(u,BoundaryCondition)
%用于求D_y^+, 默认循环边界
% uy=imfilter(u,[0;-1;1],BoundaryCondition);
%示例：uy=Dyforward(u,'circular');
if nargin<2
    BoundaryCondition='circular';
end
switch BoundaryCondition
    case 'circular'
        uy=[diff(u,1,1);u(1,:)-u(end,:)];
    case 'symmetric'
        uy=[diff(u,1,1);zeros(1,size(u,2))];
    case 'zero'
        uy=[diff(u,1,1);-u(end,:)];
end
end