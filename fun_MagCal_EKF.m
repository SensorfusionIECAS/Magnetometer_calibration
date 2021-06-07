function [dataout] = fun_MagCal_EKF(datain)
%{

���룺datain ������
    datain.h_p1 ��һʱ�̴ų���У׼ֵ
    datain.h_p0 ��ʱ�̴����ƶ���
    datain.w_cal ��ʱ�������ǵ�У׼ֵ
    datain.dt ����ʱ����
    datain.phi �����ǵĿ��Ʋ���
    datain.mrw ������������׼��
    datain.wrw ������������׼��
    datain.W ��һʱ�̵���ž���
    datain.V ��һʱ�̵�Ӳ�ž���
    datain.P ��һʱ�̵ķ������

�����dataout ������
    dataout.B ��ʱ��У׼��Ĵų�����
    dataout.W ��ʱ��У׼�����ž���
    dataout.V ��ʱ��У׼���Ӳ�ž���
    dataout.P ��ʱ��У׼��ķ������

��дʱ�䣺
  2018.4.17
%}

%% preparation
h_p_now = datain.h_p1;
h_p_last = datain.h_p0;
w_cal_last = datain.w_cal; 
parameter = datain.dt*datain.phi;
mrw = datain.mrw;
wrw = datain.wrw ;
W_last = datain.W;
b_last = datain.V;
P_last = datain.P;

%% calculate the rotarion matrix Ak
Ak = eye(3);
w_calX = [0 -w_cal_last(3) w_cal_last(2);
          w_cal_last(3) 0 -w_cal_last(1);
          -w_cal_last(2) w_cal_last(1) 0];                                  % equ 11       
Ak = Ak + w_calX * Ak * 0.01 ;                                            % equ 12b
Ak = matrix_normal(Ak);
Ak = Ak';%%

%% equation 21a
h_p_pre   = Ak * h_p_last;                                                  % equ 15a
X_pre     = [h_p_pre; W_last(1,1); W_last(2,2); W_last(3,3);                % equ 14
             W_last(1,2); W_last(1,3); W_last(2,3); b_last];                % equ 14

%% equation 21b
   Q_element = parameter*wrw.*[(h_p_last(2)+h_p_last(3)) 0 0
                            0 (h_p_last(1)+h_p_last(3)) 0
                            0 0 (h_p_last(1)+h_p_last(2))];

%    Q_element = parameter*wrw.*mean(h_p_last*h_p_last')

Qk        = [Q_element^2 zeros(3,9)
             zeros(9,3)  zeros(9,9)];                                       % equ17
Fk        = [Ak zeros(3,9);
             zeros(9,3) eye(9,9)];                                          % equ 18
P_pre     = Fk * P_last * Fk' + Qk;

yk        = h_p_now - W_last * h_p_pre - b_last;

hWk       = [h_p_pre(1) 0 0 h_p_pre(2) h_p_pre(3) 0;
             0 h_p_pre(2) 0 h_p_pre(1) 0 h_p_pre(3);
             0 0 h_p_pre(3) 0 h_p_pre(1) h_p_pre(2)];                       % equ 9b
hk        = [W_last hWk eye(3)];                                            % equ 9a

Rk        = mrw^2 * eye(3);                                                 % equ 20
Sk        = hk * P_pre * hk' + Rk;

Kk        = P_pre * hk' / Sk;

X_now     = X_pre + Kk * yk;

P_now     = (eye(12) - Kk * hk) * P_pre;

%% output the result
h_p_cal   = X_now(1:3);
W_last    = [X_now(4) X_now(7) X_now(8)
             X_now(7) X_now(5) X_now(9)
             X_now(8) X_now(9) X_now(6)];
% normalization ��һ��
balance   = abs(W_last(1,1) * W_last(2,2) * W_last(3,3))^(1/3);
% balance = trace(W_last)/3;
dataout.W = W_last / balance;
dataout.B = h_p_cal * balance;
dataout.V = X_now(10:12);
dataout.P = P_now;
dataout.h_p_pre = h_p_pre;
dataout.Kk = Kk;
ba = balance;
end
