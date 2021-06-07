function [Ae, X0, FittingResidual] = EllipsoidFitting( x, y, z, IsPlot,plottype)
% 本程序对三维散点数据进行椭球拟合
% [Ae, X0] = EllipsoidFitting( x, y, z)
% 输入参数:
%   x, y, z：分别为三维散点数据坐标;
%   IsPlot：是否作图
%   plottype: 作图类型，用于区分校正前后的球体。
% 输出参数：
%   Ae = [a d e; d b f; e f c] 为椭球的形状参数  !!并非W!!
%   X0 为椭球的中心点坐标
%   FittingResidual：拟合残差

% 二次曲面的一般形式：a x^2 + b y^2 + c z^2 + 2d xy + 2e xz + 2f yz + 2p x + 2q y + 2r z + g = 0 
% 旋转不变量：
%     I1 = a + b + c;
%     I2 = a*b + b*c + a*c -d^2 - e^2 - f^2;
%     I3 = det([a d e; d b f; e f c]);
%     I4 = det([a d e p; d b f q; e f c r; p q r g]);
% 椭球的条件：(I1 ~= 0) && (I2 > 0 ) && (I1*I3 > 0) && (I4 < 0)
% 椭球的标准形式： (X-X0)'Ae*(X-X0) = 1
%   其中：Ae = [a d e; d b f; e f c] 为椭球的形状参数
%       X0 为椭球的中心点坐标

[x0,y0,z0] = deal(x(:),y(:),z(:));
x = x(:);
y = y(:);
z = z(:);
% 数据归一化
mx = mean(x); 
my = mean(y); 
mz = mean(z); 
sx = (max(x)-min(x))/2; 
sy = (max(y)-min(y))/2; 
sz = (max(z)-min(z))/2; 
% 数据归一化
x = (x-mx)/sx; 
y = (y-my)/sy; 
z = (z-mz)/sz; 

D = [ x.^2, y.^2, z.^2, 2*x.*y, 2*x.*z, 2*y.*z, 2*x, 2*y, 2*z, ones(size(x))];
S = D'*D;

% 建立 10x10 的椭球约束矩阵 
C = diag([-1 -1 -1 zeros(1, 7)]);
I1 = [  2, 3, 11,  13,  21, 22];
I2 = [  34, 45, 56];

% 拟合成功标志，若拟合成功，Flag = 1;
Flag = 0;
k_max = 1e5;
% 判断最大的k值时，拟合二次曲面是否椭球
C(I1) = k_max/2 - 1;
C(I2) = -k_max;
[vectors,values]=eig(S, C);
values = diag(values);
I = find(real(values) > -1e-15 & ~isinf(values));
if length(I) >= 1
    [M, Index] = max(values(I));
    V = real(vectors(:,I(Index)));
    % 旋转不变量
    [a, b, c, d, e, f, p, q, r ,g] = deal(V(1), V(2), V(3), V(4), V(5), V(6), V(7) ,V(8) ,V(9) ,V(10));
    I1 = a + b + c;
    I2 = a*b + b*c + a*c -d^2 - e^2 - f^2;
    I3 = det([a d e; d b f; e f c]);
    I4 = det([a d e p; d b f q; e f c r; p q r g]);
    if (I1 ~= 0) && (I2 > 0 ) && (I1*I3 > 0) && (I4 < 0)
        % 拟合成功，获得设定限制内的最大椭球
        Flag = 1;
    end
else
    I = [];
end

if Flag == 0
    % 判断最大的k值时,拟合不成功，2分法搜索kmax
    RearchArea_min = 3;
    RearchArea_max = k_max;
    while isempty(I) && abs(RearchArea_max-RearchArea_min) > 1e-5
        k_max = (RearchArea_min + RearchArea_max)/2;
        C(I1) = k_max/2 - 1;
        C(I2) = -k_max;
        [vectors,values]=eig(S, C);
        values = diag(values);
        I = find(real(values) > -1e-15 & ~isinf(values));
        if length(I) == 1
            V = real(vectors(:,I));
            % 旋转不变量
            [a, b, c, d, e, f, p, q, r ,g] = deal(V(1), V(2), V(3), V(4), V(5), V(6), V(7) ,V(8) ,V(9) ,V(10));
            I1 = a + b + c;
            I2 = a*b + b*c + a*c -d^2 - e^2 - f^2;
            I3 = det([a d e; d b f; e f c]);
            I4 = det([a d e p; d b f q; e f c r; p q r g]);    
            if (I1 == 0) || (I2 <= 0 ) || (I1*I3 <= 0) || (I4 >= 0)
                % 拟合二次曲面不是椭球
                RearchArea_max = k_max;
                I = [];
            else
                % 拟合二次曲面是椭球
                RearchArea_min = k_max;
            end
        else
            % 拟合二次曲面不是椭球
            RearchArea_max = k_max;
            I = [];
        end
    end
end

% 反归一化
V = [V(1:3); 2*V(4:9); V(10)];
V = [   V(1)*sy^2*sz^2;
        V(2)*sx^2*sz^2;
        V(3)*sx^2*sy^2;
        V(4)*sx*sy*sz^2;
        V(5)*sx*sy^2*sz;
        V(6)*sx^2*sy*sz;
        -2*V(1)*sy^2*sz^2*mx - V(4)*sx*sy*sz^2*my - V(5)*sx*sy^2*sz*mz + V(7)*sx*sy^2*sz^2;
        -2*V(2)*sx^2*sz^2*my - V(4)*sx*sy*sz^2*mx - V(6)*sx^2*sy*sz*mz + V(8)*sx^2*sy*sz^2;
        -2*V(3)*sx^2*sy^2*mz - V(5)*sx*sy^2*sz*mx - V(6)*sx^2*sy*sz*my + V(9)*sx^2*sy^2*sz;
        V(1)*sy^2*sz^2*mx^2 + V(2)*sx^2*sz^2*my^2 + V(3)*sx^2*sy^2*mz^2 ...
          	+ V(4)*sx*sy*sz^2*mx*my + V(5)*sx*sy^2*sz*mx*mz + V(6)*sx^2*sy*sz*my*mz ...
            - V(7)*sx*sy^2*sz^2*mx - V(8)*sx^2*sy*sz^2*my - V(9)*sx^2*sy^2*sz*mz ...
            + V(10)*sx^2*sy^2*sz^2];
V = [V(1:3); V(4:9)/2; V(10)];

A = [V(1) V(4) V(5);V(4) V(2) V(6); V(5) V(6) V(3)];
b = [V(7); V(8); V(9)];
% 化为椭球的标准形式
X0 = -inv(A)*b;
Ae = A / (X0'*A*X0 - V(10));

V0 = V/(X0'*A*X0 - V(10));
D0 = [ x0.^2, y0.^2, z0.^2, 2*x0.*y0, 2*x0.*z0, 2*y0.*z0, 2*x0, 2*y0, 2*z0, ones(size(x0))];
FittingResidual = D0 * V0;

if IsPlot == 1
%     if plottype == 1
%         plot3(x0,y0,z0,'r.')
%     elseif plottype == 2
%         plot3(x0,y0,z0,'b.')
%     elseif plottype == 3
%         plot3(x0,y0,z0,'g.')
%     end
%     hold on

    % 生成标准圆球：X_sphere'*X_sphere = 1
    N_sphere = 20;
    [x_sphere, y_sphere, z_sphere] = sphere(N_sphere);
    % 椭球方程：(X-X0)'*Ae*(X-X0) = 1
    % 将Ae矩阵进行Cholesky分解，Ae = R'*R;
    % 椭球方程可写为：((X-X0)'*R')*(R *(X-X0)) = 1
    % 因此椭球点为：X_sphere = R *(X-X0);
    %   X = inv(R)*X_sphere + X0
    R = chol(Ae);
    xyz_ellipsoid = inv(R)*[x_sphere(:)'; y_sphere(:)'; z_sphere(:)'] ...
        + X0 * ones(1,(N_sphere+1)*(N_sphere+1));
    x_ellipsoid = reshape(xyz_ellipsoid(1,:), size(x_sphere));
    y_ellipsoid = reshape(xyz_ellipsoid(2,:), size(y_sphere));
    z_ellipsoid = reshape(xyz_ellipsoid(3,:), size(z_sphere));
    
    % 下边的椭球绘画，可以使椭球看起来更协调
    if plottype == 1 
        surf(x_ellipsoid, y_ellipsoid, z_ellipsoid,'FaceColor', 'r', 'EdgeColor', 'none');
    elseif plottype == 2
        surf(x_ellipsoid, y_ellipsoid, z_ellipsoid,'FaceColor', 'b', 'EdgeColor', 'none');
    elseif plottype == 3
        surf(x_ellipsoid, y_ellipsoid, z_ellipsoid,'FaceColor', 'g', 'EdgeColor', 'none');
    end
    alpha(0.4);
    axis vis3d;
    axis equal;
    camlight;
    lighting phong;
    grid on
    xlabel('X');ylabel('Y');zlabel('Z');
    
end