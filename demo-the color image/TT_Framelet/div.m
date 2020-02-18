function [u]=div(px,py,BoundaryCondition)
%用于求梯度算子的负共轭算子散度, 默认循环边界
%u=D_x^-(px)+D_y^-(py);
%示例: u=div(px,py,'circular');
if nargin<3
    BoundaryCondition='circular';
end
u=Dxbackward(px,BoundaryCondition)+Dybackward(py,BoundaryCondition);
end