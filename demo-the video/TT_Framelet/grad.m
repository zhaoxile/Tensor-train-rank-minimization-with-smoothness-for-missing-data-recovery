function [ux,uy]=grad(u,BoundaryCondition)
%用于求梯度, 默认循环边界
%ux=D_x^+(u);uy=D_y^+(u);
%示例：[ux,uy]=grad(u,'circular');
if nargin<2
    BoundaryCondition='circular';
end
ux=Dxforward(u,BoundaryCondition);
uy=Dyforward(u,BoundaryCondition);
end