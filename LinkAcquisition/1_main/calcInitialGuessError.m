%{
リンクが取れない間に伝搬する誤差を計算する
%}

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

%% setting parameter
%% constant value
constant.sunMu         = SSD.GM(10);      % km^3/s^2 gravity constant of the sun

%% parameter related to time
% simulation timeStep[s]
time.simDt = 100;
% number of time step
time.stepNum = 650; 
% simulateion start time (ephemeris time)
time.t0 = cspice_str2et('2030/01/01 00:00:00 UTC');
time.list = linspace(time.t0,time.t0+time.simDt * time.stepNum,time.stepNum+1);

%% parameter related to error
% 初期探査機軌道誤差[km]. 1000kmに設定
error.scPos0 = randn(3,1) *  1e3;
% 適当に0.1km/s程度の誤差とする
error.scVel0 = randn(3,1) *  1e-1;
% ダイナミクスの不確定性の標準偏差(探査機)
error.dynamics = 1e-10;
% 探査機状態量の初期値
sc.state0 = cspice_spkezr('699', time.t0,'ECLIPJ2000', 'NONE', '10');
      

%% 1 set parameter of s/c state
% 推定値
initialGuess = Spacecraft(time,constant.sunMu,"spacecraft");
initialGuess.state=sc.state0;

% 真値(誤差を入れる)
scTrue = Spacecraft(time,constant.sunMu,"spacecraft");
scTrue.state = initialGuess.state + [error.scPos0; error.scVel0];

%% 2. calculate orbit of spacecraft
scTrue.calcOrbitTwoBody(time.simDt,error.dynamics)
initialGuess.calcOrbitTwoBody(time.simDt,0)



% 軌道決定精度をplotする
% 探査機位置誤差の時間履歴
figure(1)
title('position error of estimated value')
hold on
plot((time.list-time.list(1))/60/60/24, initialGuess.state(1,:) - scTrue.state(1,:))
plot((time.list-time.list(1))/60/60/24, initialGuess.state(2,:) - scTrue.state(2,:))
plot((time.list-time.list(1))/60/60/24, initialGuess.state(3,:) - scTrue.state(3,:))
plot((time.list-time.list(1))/60/60/24,...
    ( (initialGuess.state(1,:) - scTrue.state(1,:)).^2 + (initialGuess.state(2,:) - scTrue.state(2,:)).^2 + (initialGuess.state(3,:) - scTrue.state(3,:)).^2).^0.5)
xlabel('day')
ylabel('position error [km]')
legend('x', 'y', 'z','length')
hold off

figure(2)
title('velocity error of estimated value')
hold on
plot((time.list-time.list(1))/60/60/24, initialGuess.state(4,:) - scTrue.state(4,:))
plot((time.list-time.list(1))/60/60/24, initialGuess.state(5,:) - scTrue.state(5,:))
plot((time.list-time.list(1))/60/60/24, initialGuess.state(6,:) - scTrue.state(6,:))
plot((time.list-time.list(1))/60/60/24,...
    ( (initialGuess.state(4,:) - scTrue.state(4,:)).^2 + (initialGuess.state(5,:) - scTrue.state(5,:)).^2 + (initialGuess.state(6,:) - scTrue.state(6,:)).^2).^0.5)
xlabel('day')
ylabel('position error [km]')
legend('x', 'y', 'z','length')
hold off
