# Magnetometer_calibration
use EKF and IMU(gyroscope)  to calibration the magnetometer

# 主程序
ResultCompare
可直接运行

1. EllipsoidFitting.m 本程序对三维散点数据进行椭球拟合
2. f_Mag_Ellipsoid_Fit.m 椭球拟合函数。包括四分量、七分量和十分量。
3. fun_MagCal_EKF.m EKF方法
4. huawei_x1.xlsx 使用华为手机采集的数据集
