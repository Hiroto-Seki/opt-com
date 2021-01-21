%% main関数 optComAcquisition_main
%% 前処理
clear all; close all; clc
% add path to SPICE
addpath(genpath('~/Documents/Matlab/SPICE'));
% SPICEのKernel(天体情報)を読み込む
spice_loadkernels();
% 取得した天体情報 + alphaを利用しやすいように構造体へまとめる
SSD = spice_setparams();
% 乱数
rng('default');
rng(2)

[constant,time,error,gs,sc,gsTrue,earth,scTrue,scEstByScEkf,scEstByGsEkf,ekf,~] = setparam(SSD);
% 長時間の通信隔絶前の捕捉
[gsTrue,earth,scTrue,scEstByScEkf,scEstByGsEkf,time] = optComAcquisition_EKF(constant,time,error,gs,sc,gsTrue,earth,scTrue,scEstByScEkf,scEstByGsEkf,ekf,SSD);
showResult(scTrue,scEstByScEkf,scEstByGsEkf,error,0);
% 長時間の通信隔絶
timeUpdateDuringLongLos;
% 長時間の通信隔絶後の再捕捉
[time_afterLos,gsTrue_afterLos,earth_afterLos,scTrue_afterLos,scEstByScEkf_afterLos,scEstByGsEkf_afterLos]...
    = setParam_afterLos(time,time_los,gs,scTrue_los,scEstBySc_los,scEstByGs_los,scEstByScEkf,scEstByGsEkf,constant,error);
[gsTrue_afterLos,earth_afterLos,scTrue_afterLos,scEstByScEkf_afterLos,scEstByGsEkf_afterLos,time_afterLOS] = optComAcquisition_EKF(constant,time_afterLos,error,gs,sc,gsTrue_afterLos,earth_afterLos,scTrue_afterLos,scEstByScEkf_afterLos,scEstByGsEkf_afterLos,ekf,SSD);
showResult(scTrue_afterLos,scEstByScEkf_afterLos,scEstByGsEkf_afterLos,error,2);
