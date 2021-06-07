%{
分别知识EKF校正和椭球拟合，比较二者的效果。

编写时间：
  2018.4.18
%}
%% format
clear

%% load the data
% [gyro_still] = xlsread('./stableB.xlsx');                                   % 静止时陀螺仪的输出！！！
% figure
% plot(gyro_still(1:1000,4:6))
[data_test] = xlsread('huawei_x1.xlsx');
mag_test = data_test(:,1:3);

% [data_test]  = xlsread ('数据集\fastwalking_swing_circle.xlsx');     
% mag_test = data_test(:,1:3)/10;
testlength      = size( data_test, 1 );

% [M]          = xlsread ('数据集\fastwalking_swing_circle.xlsx');                          % 载入初始数据
% data.b_p     = M(:,1:3)/10;                                                    % 磁力计测量值，单位为mG
% data.w       = M(:,4:6)/1800*pi;           
[M]          = xlsread ('huawei_x1.xlsx');                          % 载入初始数据
data.b_p     = M(:,1:3);                                                    % 磁力计测量值，单位为mG
data.w       = M(:,4:6);                                                    % 陀螺仪测量值，单位rad/s
data.dt      = 0.01;                                                        % 量测时间间隔，单位s
data.mrw     = 0.5;                                                         % 磁场b的random walks，单位为mG，依照论文中table I设置的值
data.wrw     = 0.1/180*pi;                                                  % 角速度w的random walks,单位为degree/（s^1/2）依照论文中table I设置的值                                     % 测量次数
data.m       = size( data.b_p, 1 );
data.phi     = 100;                                                         % 论文中公式17中的参数，用于表示陀螺仪的可靠程度。
data.P       = [500*eye(3) zeros(3,6)  zeros(3);
                zeros(6,3) 1e-4*eye(6) zeros(6,3);
                zeros(3)   zeros(3,6)  500*eye(3)];                         % 方差矩阵，依照论文中table I设置的值
W = eye(3);
V = [0 0 0]';
%% calibrated gyroscope measurement at time k-1，equ.16
% w_bias     = mean(gyro_still(1:1000,4:6))/1800*pi;
% for j = 1: data.m
%     w_noise(j,:)     = w_bias; 
% end
w            = data.w;
sz_w         = size(w);
w_noise      = data.wrw * randn(sz_w); 
w_cal        = data.w - w_noise;

%% Select correction interval
% 5s movement data for calibration
start       = 500;
range       = start+1:1:start+1+1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EKF Calibration
k = 1;
for i = range 
    %% make the datain
    datain.h_p1       = data.b_p(i+1,:)';
    datain.h_p0       = data.b_p(i,:)';
    datain.dt         = data.dt;
    datain.phi        = data.phi;
    datain.mrw        = data.mrw;
    datain.wrw        = data.wrw;   
    datain.w_cal      = w_cal(i,:);
    if k == 1 
    % Setting the initial value
        datain.W      = eye(3);
        datain.V      = [0 0 0]';
        datain.P      = data.P;
        datain.h_p0   = datain.W ^(-1) * (datain.h_p0 - datain.V);
    else
    % make the value equal to the value at last time 
        datain.W      = dataout.W;
        datain.V      = dataout.V;
        datain.P      = dataout.P;
        datain.h_p0   = dataout.B;
    end
    
    %% EKF
    [dataout] = fun_MagCal_EKF(datain);
    
    %% load the output
    B_cal_EKF(k,:)    = dataout.B;
    W_cal_EKF(:,:,k)  = [dataout.W(1,1) dataout.W(1,2) dataout.W(1,3);
                         dataout.W(1,2) dataout.W(2,2) dataout.W(2,3);
                         dataout.W(1,3) dataout.W(2,3) dataout.W(3,3)];
    V_cal_EKF(k,:)    = dataout.V;
    P_cal_EKF(k,:)    = [dataout.P(1,1) dataout.P(4,4) dataout.P(5,5) dataout.P(10,10)];
    h_pre_EKF(k,:)    = dataout.h_p_pre;
    Kk(:,:,k)    = dataout.Kk;

    k = k + 1;
end
%% output
EKF_W = W_cal_EKF(:,:,end);
EKF_V =V_cal_EKF(end,:)';

figure 
plot(V_cal_EKF(:,1))
% 显示结果
fprintf( '软磁矩阵 EKF_W:\n  %g %g %g\n  %g %g %g\n  %g %g %g\n',EKF_W);
fprintf( '硬磁矢量 EKF_V:\n  %g %g %g\n',EKF_V);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ellipsoid Fitting
h = 1;
for i = range 
    dataElli(h,:) = data.b_p(i,:);
    h = h + 1; 
end
[Elli_Winv, Elli_V, Elli_B, Elli_E] = f_Mag_Ellipsoid_Fit(dataElli,10);
% 显示结果
fprintf( '软磁矩阵 Elli_Winv:\n  %g %g %g\n  %g %g %g\n  %g %g %g\n',Elli_Winv);
fprintf( '硬磁矢量 Elli_V:\n  %g %g %g\n',Elli_V);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compare
for i= 1:testlength
    normdata(i,:)      = norm(mag_test(i,:));

    B_EKF_WV(i,:)      = EKF_W^(-1)*(mag_test(i,:)' - EKF_V);
    normdata_EKF(i,:)  = norm(B_EKF_WV(i,:));

    B_Elli_WV(i,:)     = Elli_Winv * (mag_test(i,:)' - Elli_V);
    normdata_Elli(i,:) = norm(B_Elli_WV(i,:));

end


% figure
% plot3(h_pre_EKF(:,1),h_pre_EKF(:,2),h_pre_EKF(:,3),'.')
% hold on
% plot3(data.b_p(range,1),data.b_p(range,2),data.b_p(range,3),'.')
% hold on 
% % plot3(B_cal_EKF(:,1),B_cal_EKF(:,2),B_cal_EKF(:,3),'.')
% axis equal
% grid on 



%% comppute the magnetic field strength
Bfield_EKF  = mean(normdata_EKF);
Bfield_Elli = mean(normdata_Elli);
Bfield_raw      = mean(normdata);
standard = 50 * ones(testlength ,1);

quality1_EKF  = std(normdata_EKF) / Bfield_EKF
quality1_Elli = std(normdata_Elli) / Bfield_Elli

difference_raw = normdata-standard;
difference_EKF = normdata_EKF-standard;
difference_Elli = normdata_Elli-standard;
figure
plot(difference_raw,'r');hold on 
plot(difference_EKF,'b');hold on 
plot(difference_Elli,'g');hold on 
legend('原始数据','EKF后数据','Elli后数据')
xlabel('时间(unit:0.01s)')
ylabel('距离(unit:mG)')



quality2_raw  = mean(normdata-standard)
quality2_EKF  = mean(normdata_EKF-standard)
quality2_Elli = mean(normdata_Elli-standard)

figure
grid on
% plot3(data.b_p(:,1),data.b_p(:,2),data.b_p(:,3),'.')
% hold on
 plot3(B_EKF_WV(:,1),B_EKF_WV(:,2),B_EKF_WV(:,3),'.')
% hold on
[x, y, z] = ellipsoid(0,0,0,Bfield_EKF,Bfield_EKF,Bfield_EKF);
mesh(x, y, z, 'FaceAlpha', 0.1, 'EdgeColor', 'k');
hold on
[~, X0, ~] = EllipsoidFitting(mag_test(:,1),mag_test(:,2),mag_test(:,3), 1, 1);
[~, X_EKF, ~] = EllipsoidFitting( B_EKF_WV(:,1) , B_EKF_WV(:,2) , B_EKF_WV(:,3), 1, 2);
[~, X_Elli, ~] = EllipsoidFitting( B_Elli_WV(:,1),B_Elli_WV(:,2),B_Elli_WV(:,3), 1, 3);
legend('理想磁场球体','原始椭球','EKF椭球','Elli椭球');
% legend('理想磁场球体','原始椭球','EKF椭球');

xlabel('X (unit:mG)')
ylabel('Y (unit:mG)')
zlabel('Z (unit:mG)')
plot3(data.b_p(:,1),data.b_p(:,2),data.b_p(:,3),'r.');hold on
plot3(data.b_p(range,1),data.b_p(range,2),data.b_p(range,3),'ys');hold on
axis vis3d;
axis equal;
grid on
% title('效果展示')
% axis([-80 80 -80 80 -80 80 ])

% axis equal;
% grid on
% title('Elli')
% [x_Elli, y_Elli, z_Elli] = ellip% figure
% plot3(data.b_p(:,1),data.b_p(:,2),data.b_p(:,3),'.')
% hold on
% plot3(B_Elli_WV(:,1),B_Elli_WV(:,2),B_Elli_WV(:,3),'.')
% axis vis3d;
% ellipsoid(0,0,0,Bfield_Elli,Bfield_Elli,Bfield_Elli);
% mesh(x_Elli, y_Elli, z_Elli, 'FaceAlpha', 0.1, 'EdgeColor', 'r');





