% 複数回のシミュレーション結果をまとめる

clear all; close all; clc
rng('default');

% シミュレーション回数
simNum = 10;

% 結果のリスト

% 初期時刻での探査機自身の位置推定値の誤差(スカラー)
resultList.scEstPosError0 = zeros(1,length(simNum));
% 最終時刻での探査機自身の位置推定値の誤差(スカラー)
resultList.scEstPosError = zeros(1,length(simNum));

% 初期時刻での探査機自身の速度推定値の誤差(スカラー)
resultList.scEstVelError0 = zeros(1,length(simNum));
% 最終時刻での探査機自身の速度推定値の誤差(スカラー)
resultList.scEstVelError = zeros(1,length(simNum));

% 最終時刻での探査機によるクロック誤差
resultList.scEstClockError0 =  zeros(1,length(simNum));
% 最終時刻での探査機によるクロック誤差
resultList.scEstClockError =  zeros(1,length(simNum));

% 初期時刻での地上局による探査機の位置推定値の誤差(スカラー)
resultList.gsEstPosError0 = zeros(1,length(simNum));
% 最終時刻での地上局による探査機の位置推定値の誤差(スカラー)
resultList.gsEstPosError = zeros(1,length(simNum));

% 初期時刻での地上局による探査機の速度推定値の誤差(スカラー)
resultList.gsEstVelError0 = zeros(1,length(simNum));
% 最終時刻での地上局による探査機の速度推定値の誤差(スカラー)
resultList.gsEstVelError = zeros(1,length(simNum));


for h = 1:simNum
    rng(h)
    LinkAcquisition();
    LastID = length(time.list);
    % 結果を格納
    resultList.scEstPosError0(h)   = norm(scEst.state(1:3,1) - scTrue.state(1:3,1));
    resultList.scEstVelError0(h)   = norm(scEst.state(4:6,1) - scTrue.state(4:6,1));
    resultList.scEstClockError0(h) = scEst.resClockError(1);
    resultList.gsEstPosError0(h)   = norm(scEstGs.state(1:3,1) - scTrue.state(1:3,1));
    resultList.gsEstVelError0(h)   = norm(scEstGs.state(4:6,1) - scTrue.state(4:6,1));
     
    resultList.scEstPosError(h)   = norm(scEst.state(1:3,LastID) - scTrue.state(1:3,LastID));
    resultList.scEstVelError(h)   = norm(scEst.state(4:6,LastID) - scTrue.state(4:6,LastID));
    resultList.scEstClockError(h) = scEst.resClockError(LastID);
    resultList.gsEstPosError(h)   = norm(scEstGs.state(1:3,LastID) - scTrue.state(1:3,LastID));
    resultList.gsEstVelError(h)   = norm(scEstGs.state(4:6,LastID) - scTrue.state(4:6,LastID));
    
end


% 初期誤差の平均値
result.PosError0 = mean(resultList.scEstPosError0);
result.VelError0 = mean(resultList.scEstVelError0);

