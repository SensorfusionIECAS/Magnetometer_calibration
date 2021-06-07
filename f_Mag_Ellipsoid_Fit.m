function [ W, V, B, E ] = f_Mag_Ellipsoid_Fit( X, flag)
%{
������Ϻ����������ķ������߷�����ʮ������
���룺
  X, [x y z]    - У�����ݣ�Ϊn * 3�ľ���
  flag          - 4 ���ķ����������
                - 7 ���߷����������
                - 10 ��ʮ�����������
�����
  W             - ���У������
  V             - Ӳ��У��ʸ��
  B             - �شų�ǿ��
  E             - ������ ��e = 0.5*sqrt(E/M),��������ο����µڶ������ݣ�
�ο���
  ��Magnetic Calibration��
��дʱ�䣺
  2018.2.8
%}
%% �������
narginchk( 1, 2); % ������뺯���Ĳ�������
  % �������ֵ�ø�������1-2֮�䣬�������������ڣ����κη�Ӧ��
  
if nargin == 1 % ��������������Ϊ1��Ĭ��ʹ��ʮ�����������
    flag = 10;  
end

if size( X, 2 ) ~= 3 % �������Ĵ���������
    error( 'Input data must have three columns!' );
else
    x = X( :, 1 );
    y = X( :, 2 );
    z = X( :, 3 );
end

%% ��ʼ������
W = eye(3);  
V = [1 1 1];
B = 0;
E = 0;
m = size( x, 1 ); % ��������

%% �ķ������
    % Normal Equations Solution in Non-Homogeneous Case���������µ��������̽�
if flag == 4
    XX = [ x, y, z, ones( m, 1 )];
    YY = x.^2 + y.^2 + z.^2;
    s = ( XX' * XX )\( XX' * YY );
    
    V = (1/2) * [s(1) s(2) s(3)]';
    B = sqrt(s(4) + V(1)^2 + V(2)^2 + V(3)^2);
    % E - ���㸴�ӣ���ʱû�б�д
    % W - WΪ��λ��

%% �߷������
elseif flag == 7
    % Eigenvector Solution in Homogeneous Case�������µ�����ֵ�����
    XX = [ x.^2, y.^2, z.^2, x, y, z, ones( m, 1 )];
    [vector ,value_matrix] = eig( XX' * XX);
    value = value_matrix(  value_matrix ~= 0 );
    minposition = find( value == min( value ));
    s = vector( :,minposition );
    A = [ s(1) 0 0 ; 0 s(2) 0 ; 0 0 s(3)];
    
    % ���A������ʽ�Ƿ�С���㣬���С���㣬ʹ��s���
    d = det(A); % ����A������ʽ
    if d < 0
        s = -s;
        A = -A; % ����A��ֵ
        d = -d;
    end
    
    % ����A������ʽ���Խ�s�������š���������ֱ˵�˶�A�������ţ�������ܲ��ã�Ӧ�ò��ԣ�    
    scalar = d^(-1/3);
    s = s * scalar;
    
    % ����A��ֵ
    A = [ s(1) 0 0 ; 0 s(2) 0 ; 0 0 s(3)];

    % ������
    B =  sqrt( abs( A(1,1)*V(1)^2 + A(2,2)*V(2)^2 + A(3,3)*V(3)^2 - s(7)));
    V = [ -s(4)/(2 * s(1)), -s(5)/(2 * s(2)), -s(6)/(2 * s(3))]';
    W = sqrt( abs( A));
    E = ( 1 / ( 2*B^2)) * sqrt( value( minposition) / m);
    
%% ʮ�������
elseif flag == 10
    XX = [ x.^2, 2.*x.*y, 2.*x.*z, y.^2, 2.*y.*z, z.^2, x, y, z, ones( m ,1)];
    [vector ,value_matrix] = eig( XX' * XX);
    value = value_matrix(  value_matrix ~= 0 );
    minposition = find( value == min( value ));
    s = vector( :,minposition );
    
    % ���A������ʽ�Ƿ�С���㣬���С���㣬ʹ��s���
    A = [ s(1) s(2) s(3) ;...
          s(2) s(4) s(5) ;...
          s(3) s(5) s(6)]; 
    d = det(A);
    if d < 0
        s = -s;
        d = -d;
    end
    
    % ����A������ʽ���Խ�s�������š���������ֻ˵�˶�A�������ţ�������ܲ��ã�Ӧ�ò��ԣ�
    scalar= d^(-1/3);
    s = s * scalar;
    
    % ����A��ֵ
    A = [ s(1) s(2) s(3) ;...
          s(2) s(4) s(5) ;...
          s(3) s(5) s(6)]; 
      
    % ������
    V = (-0.5) * (A \ [s(7) s(8) s(9)]'); % A�е����ź�S��7.8.9���е����������������V����
    B = sqrt( abs( A(1,1)*V(1)^2 + 2*A(1,2)*V(1)*V(2) + 2*A(1,3)*V(1)*V(3)...
                 + A(2,2)*V(2)^2 + 2*A(2,3)*V(2)*V(3) + A(3,3)*V(3)^2 - s(10)));
    [vectorA ,valueA] = eig( A);
    W = vectorA * sqrt( abs( valueA)) * vectorA';
    % W�ļ��������������ַ������˴���������ֵ������Ϊ����ֱ�ӿ����ţ��õ��Ľ�����á�
    E = ( 1 / ( 2*B^2)) * sqrt( value( minposition) / m);
end

%% ��ʾ���
% fprintf( '���� %g ����������Ͻ��дų�У�������Ϊ:\n', flag );
fprintf( 'Ӳ��У��ʸ�� V:\n  %g %g %g\n', V );
% fprintf( '���У������ W:\n  %g %g %g\n  %g %g %g\n  %g %g %g\n', W );
fprintf( '�شų�ǿ�� B:\n  %g \n', B );
% fprintf( '������ E:\n  %g \n', E );

end














