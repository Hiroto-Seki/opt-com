%% main関数 optComAcquisition_main
%% 前処理
clear all; close all; clc


% 環境によってパスを変更する

if ismac
    % add path to SPICE
    addpath(genpath('~/Documents/Matlab/SPICE'));
    % add path to keep the result
    resultFile = char(datetime('now','Format','yyyyMMddHH'));
    mkdir('~/Documents/lab/master/research/result', resultFile)
    resultPath = ['~/Documents/lab/master/research/result/',resultFile];
    addpath(resultPath);
elseif ispc
    % add path to SPICE
    addpath(genpath('~/Documents/matlab/SPICE'));
    % add path to keep the result
    resultFile = char(datetime('now','Format','yyyyMMddHH'));
    mkdir('~/Documents/opt-com/result', resultFile)
    resultPath = ['~/Documents/opt-com/result/',resultFile];
    addpath(resultPath);   
end

% SPICEのKernel(天体情報)を読み込む
spice_loadkernels();
% 取得した天体情報 + alphaを利用しやすいように構造体へまとめる
SSD = spice_setparams();
% 乱数
rng('default');

% 結果の格納
N = 1; %モンテカルロシミュレーションの数

% result(1):通信断絶前 result(2):通信断絶中 K=3:通信断絶後
result = struct("posErrorSc"  , zeros(8,N),... %x,y,z,rの誤差とそれに対応する標準偏差 
                "velErrorSc"  , zeros(8,N),... %x,y,z,rの誤差とそれに対応する標準偏差
                "clockErrorSc",zeros(2,N),...  %クロックの推定誤差とそれに対応する標準偏差
                "posErrorGs"  , zeros(8,N),... %x,y,z,rの誤差とそれに対応する標準偏差 
                "velErrorGs"  , zeros(8,N),... %x,y,z,rの誤差とそれに対応する標準偏差
                "clockErrorGs",zeros(2,N),...  %クロックの推定誤差とそれに対応する標準偏差 
                "downAvail" ,zeros(1,N));    %Downlinkの成功率


for n = 3
    rng(n)

    [constant,time,error,gs,sc,gsTrue,earth,scTrue,scEstByScEkf,scEstByGsEkf,ekf,~] = setparam(SSD);
    % 長時間の通信隔絶前の捕捉
    [gsTrue,earth,scTrue,scEstByScEkf,scEstByGsEkf,time] = optComAcquisition_EKF(constant,time,error,gs,sc,gsTrue,earth,scTrue,scEstByScEkf,scEstByGsEkf,ekf,SSD);
    showResult(scTrue,scEstByScEkf,scEstByGsEkf,error,n*3-3,gsTrue,gs,resultPath);

    % 長時間の通信隔絶
    timeUpdateDuringLongLos;
    % 長時間の通信隔絶後の再捕捉
    [time_afterLos,gsTrue_afterLos,earth_afterLos,scTrue_afterLos,scEstByScEkf_afterLos,scEstByGsEkf_afterLos]...
        = setParam_afterLos(time,time_los,gs,scTrue_los,scEstBySc_los,scEstByGs_los,scEstByScEkf,scEstByGsEkf,constant,error);
    [gsTrue_afterLos,earth_afterLos,scTrue_afterLos,scEstByScEkf_afterLos,scEstByGsEkf_afterLos,time_afterLOS] = optComAcquisition_EKF(constant,time_afterLos,error,gs,sc,gsTrue_afterLos,earth_afterLos,scTrue_afterLos,scEstByScEkf_afterLos,scEstByGsEkf_afterLos,ekf,SSD);
    showResult(scTrue_afterLos,scEstByScEkf_afterLos,scEstByGsEkf_afterLos,error,n*3-1,gsTrue_afterLos,gs,resultPath);
    
    % 結果の格納
    saveResult;
end

save([resultPath,'/result.mat'],'result')
% comparison;






