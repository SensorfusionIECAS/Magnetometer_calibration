function [ W, V, B, E ] = f_Mag_Ellipsoid_Fit( X, flag)
%{
椭球拟合函数。包括四分量、七分量和十分量。
输入：
  X, [x y z]    - 校正数据，为n * 3的矩阵。
  flag          - 4 ，四分量椭球拟合
                - 7 ，七分量椭球拟合
                - 10 ，十分量椭球拟合
输出：
  W             - 软磁校正矩阵
  V             - 硬磁校正矢量
  B             - 地磁场强度
  E             - 拟合误差 （e = 0.5*sqrt(E/M),具体意义参考文章第二章内容）
参考：
  “Magnetic Calibration”
编写时间：
  2018.2.8
%}
%% 检查输入
narginchk( 1, 2); % 检查输入函数的参数个数
  % 如果输入值得个数不在1-2之间，则输出错误，如果在，无任何反应。
  
if nargin == 1 % 如果输入变量个数为1，默认使用十分量椭球拟合
    flag = 10;  
end

if size( X, 2 ) ~= 3 % 检查输入的磁力计数据
    error( 'Input data must have three columns!' );
else
    x = X( :, 1 );
    y = X( :, 2 );
    z = X( :, 3 );
end

%% 初始化变量
W = eye(3);  
V = [1 1 1];
B = 0;
E = 0;
m = size( x, 1 ); % 测量次数

%% 四分量拟合
    % Normal Equations Solution in Non-Homogeneous Case非齐次情况下的正常方程解
if flag == 4
    XX = [ x, y, z, ones( m, 1 )];
    YY = x.^2 + y.^2 + z.^2;
    s = ( XX' * XX )\( XX' * YY );
    
    V = (1/2) * [s(1) s(2) s(3)]';
    B = sqrt(s(4) + V(1)^2 + V(2)^2 + V(3)^2);
    % E - 计算复杂，暂时没有编写
    % W - W为单位阵。

%% 七分量拟合
elseif flag == 7
    % Eigenvector Solution in Homogeneous Case齐次情况下的特征值法求解
    XX = [ x.^2, y.^2, z.^2, x, y, z, ones( m, 1 )];
    [vector ,value_matrix] = eig( XX' * XX);
    value = value_matrix(  value_matrix ~= 0 );
    minposition = find( value == min( value ));
    s = vector( :,minposition );
    A = [ s(1) 0 0 ; 0 s(2) 0 ; 0 0 s(3)];
    
    % 检查A的行列式是否小于零，如果小于零，使解s变号
    d = det(A); % 计算A的行列式
    if d < 0
        s = -s;
        A = -A; % 更新A的值
        d = -d;
    end
    
    % 根据A的行列式，对解s进行缩放。（文章中直说了对A进行缩放，但结果很不好，应该不对）    
    scalar = d^(-1/3);
    s = s * scalar;
    
    % 更新A的值
    A = [ s(1) 0 0 ; 0 s(2) 0 ; 0 0 s(3)];

    % 输出结果
    B =  sqrt( abs( A(1,1)*V(1)^2 + A(2,2)*V(2)^2 + A(3,3)*V(3)^2 - s(7)));
    V = [ -s(4)/(2 * s(1)), -s(5)/(2 * s(2)), -s(6)/(2 * s(3))]';
    W = sqrt( abs( A));
    E = ( 1 / ( 2*B^2)) * sqrt( value( minposition) / m);
    
%% 十分量拟合
elseif flag == 10
    XX = [ x.^2, 2.*x.*y, 2.*x.*z, y.^2, 2.*y.*z, z.^2, x, y, z, ones( m ,1)];
    [vector ,value_matrix] = eig( XX' * XX);
    value = value_matrix(  value_matrix ~= 0 );
    minposition = find( value == min( value ));
    s = vector( :,minposition );
    
    % 检查A的行列式是否小于零，如果小于零，使解s变号
    A = [ s(1) s(2) s(3) ;...
          s(2) s(4) s(5) ;...
          s(3) s(5) s(6)]; 
    d = det(A);
    if d < 0
        s = -s;
        d = -d;
    end
    
    % 根据A的行列式，对解s进行缩放。（文章中只说了对A进行缩放，但结果很不好，应该不对）
    scalar= d^(-1/3);
    s = s * scalar;
    
    % 更新A的值
    A = [ s(1) s(2) s(3) ;...
          s(2) s(4) s(5) ;...
          s(3) s(5) s(6)]; 
      
    % 输出结果
    V = (-0.5) * (A \ [s(7) s(8) s(9)]'); % A中的缩放和S（7.8.9）中的缩放想抵消，所以V不变
    B = sqrt( abs( A(1,1)*V(1)^2 + 2*A(1,2)*V(1)*V(2) + 2*A(1,3)*V(1)*V(3)...
                 + A(2,2)*V(2)^2 + 2*A(2,3)*V(2)*V(3) + A(3,3)*V(3)^2 - s(10)));
    [vectorA ,valueA] = eig( A);
    W = vectorA * sqrt( abs( valueA)) * vectorA';
    % W的计算在文中有两种方法，此处采用特征值法，因为发现直接开根号，得到的结果不好。
    E = ( 1 / ( 2*B^2)) * sqrt( value( minposition) / m);
end

%% 显示结果
% fprintf( '采用 %g 分量椭球拟合进行磁场校正，结果为:\n', flag );
fprintf( '硬磁校正矢量 V:\n  %g %g %g\n', V );
% fprintf( '软磁校正矩阵 W:\n  %g %g %g\n  %g %g %g\n  %g %g %g\n', W );
fprintf( '地磁场强度 B:\n  %g \n', B );
% fprintf( '拟合误差 E:\n  %g \n', E );

end














