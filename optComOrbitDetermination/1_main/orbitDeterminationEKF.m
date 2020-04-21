% ---  Contents   ---
% orbit determination using EKF
% ---  outline    ---
% 1. get ephemris of gorund station 
% 2. calculate true orbit 
% 3. calculate observed value using true value
% 4. estimate spacecraft orbit using observed value by EKF
% ---   Date      ---
% rev1: 2020/04/19

clear all; close all; clc
%% setting parameter
rng('default');
rng(1)
% add path to SPICE
addpath(genpath('~/Documents/Matlab/SPICE'));
% SPICEのKernel(天体情報)を読み込む
spice_loadkernels();
% 取得した天体情報 + alphaを利用しやすいように構造体へまとめる
SSD = spice_setparams();

%% constant value
constant.sunMu         = SSD.GM(10);      % km^3/s^2 gravity constant of the sun
constant.earthMu       = SSD.GM(399);     % km^3/s^2 gravity constant of the 
constant.saturnMu      = SSD.GM(799);     %[km/s^2]
constant.lightSpeed    = 299792.458;      % [km/s] light speed 
constant.earthRadius   = SSD.r(399);      % [km] earth radius
constant.earthRotation = 2*pi/24/60/60;   % [rad/s] earth rotation speed
constant.eathAxis      = 23.4/180*pi;     % inclination of the earth axis

%% parameter related to time
% simulation timeStep[s]
time.simDt = 10;
% number of time step
time.stepNum = 1200; 
% simulateion start time (ephemeris time)
time.t0 = cspice_str2et('2030/01/01 00:00:00 UTC');
time.list = linspace(time.t0,time.t0+time.simDt * time.stepNum,time.stepNum+1);

%% 誤差に関するパラメーター
% 初期時計誤差
error.clock = 10   * randn; %時計誤差(秒)
% 初期探査機軌道誤差[km]. 地球からみて1μradくらいの誤差にする
error.scPos0 = randn(3,1) *  150;
% 適当に0.1km/s程度の誤差とする
error.scVel0 = randn(3,1) *  0.1;
% ダイナミクスの不確定性の標準偏差(探査機)
error.dynamics = 1e-10;
% 角度観測の不確定性の標準偏差(もしかしたらもっと悪いかも)
error.angle = 1e-6;

%% set gorund station (USUDAにする)
% 緯度経度から，ECLIPJ2000,地球中心座標系での座標を得たい
gs.lat  = 36.1325063*cspice_rpd();
gs.lon =138.3607113*cspice_rpd();
gs.alt  = 1.456;
constant.earthF = 1/298.257223560; % https://topex.ucsd.edu/geodynamics/14gravity1_2.pdf
gs.bodyFixed = cspice_georec(gs.lat,gs.lon,gs.alt,constant.earthRadius,constant.earthF);
gs.bodyFixed = [gs.bodyFixed;0;0;0];
ephemData.xform_ec2IAU = cspice_sxform('IAU_EARTH','ECLIPJ2000',time.list);

%% set spacecraft initial state (土星と同じ位置とする)
scEst = orbitalState('spacecraft','estimate','SolarSystem');
scEst.state0=cspice_spkezr('699', time.t0,'ECLIPJ2000', 'NONE', '10');
scEst.pos0 =  scEst.state0(1:3);
scEst.vel0 =  scEst.state0(4:6);
scTrue = orbitalState('spacecraft','true','SolarSystem');
scTrue.pos0 = scEst.pos0 + error.scPos0;
scTrue.vel0 = scEst.vel0 + error.scVel0;

%% 1. get ephemris data of ground station
% get ephemeris data of the earth orbit in ECLIPJ2000, sun center
ephemData.earth = cspice_spkezr('399', time.list,'ECLIPJ2000', 'NONE', '10');
% get ehemeris data of the ground station in ECLIPJ2000, earth center
ephemData.gs = zeros(6,time.stepNum+1);
for i_1 = 1: time.stepNum+1
    ephemData.gs(:,i_1) = ephemData.xform_ec2IAU(:,:,i_1)*gs.bodyFixed;
end

%% 2. calculate true orbit of spacecraft
scTrue.calcOrbitTwoBody_pertubation(time,constant, error.dynamics)
% 観測もなかった時の軌道も作っておく(初期値は推定値と同じ)
scInitialGuess = orbitalState('spacecraft','estimate','SolarSystem');
scInitialGuess.pos0 = scEst.pos0;
scInitialGuess.vel0 = scEst.vel0;
scInitialGuess.calcOrbitTwoBody_pertubation(time,constant, 0)

%% 3. calculate observed value, 4. estimate spacecraft orbit using observed value by EKF
%% 観測する(観測は真値を伝搬して求めていく)
observe = observeState(time);
% 推定値の初期化
scEst.clockErrorCorrection= zeros(1,time.stepNum+1);
scEst.pos = zeros(3,time.stepNum+1);
scEst.vel = zeros(3,time.stepNum+1);
% 初期推定値を入力
% 初期予想時計誤差は0にする. X_hat初期値の第一項
X_hat = [0;scEst.pos0;scEst.vel0];
   scEst.clockErrorCorrection(1)= X_hat(1);
   scEst.pos(:,1) = X_hat(2:4);
   scEst.vel(:,1) = X_hat(5:7);
P = 1e4 * eye(7);  %正直適当
R = [100,0,0;      %観測誤差
    0,1e-12,0;
    0,0,1e-12];

for i_2 = 2:length(time.list)
   %% まずは観測量を得る
   observe.calcObservedValue(time,ephemData,scTrue,constant,error,i_2)
   Y = [observe.ltd(i_2); observe.azimuth(i_2); observe.elevation(i_2)];
%    % 観測量からレーザーが発された時刻の地上局の状態量を得る
%    xve = observe.xve(:,i_2);
%    xvg = observe.xvg(:,i_2);
   %% リファレンスの状態量とSTMを伝搬して求める(RK4,1stepで積分した)
   STM = eye(7);
   xvsc = X_hat(2:7);
   k1sc = orbitalState.twobody(xvsc,constant.sunMu,0);
   k2sc = orbitalState.twobody(xvsc+0.5*time.simDt*k1sc,constant.sunMu,0);
   k3sc = orbitalState.twobody(xvsc+0.5*time.simDt*k2sc,constant.sunMu,0);
   k4sc = orbitalState.twobody(xvsc+time.simDt*k3sc,constant.sunMu,0);
%    xvsc = xvsc + time.simDt/6*(k1sc+2*k2sc+2*k3sc+k4sc);
%%   STMも時間の関数なので,RK4内でも時間の関数にしたい
   k1stm = orbitalState.delFdelX(xvsc,constant.sunMu)*STM;
   k2stm = orbitalState.delFdelX(xvsc+0.5*time.simDt*k1sc,constant.sunMu)*(STM+0.5*time.simDt*k1stm);
   k3stm = orbitalState.delFdelX(xvsc+0.5*time.simDt*k2sc,constant.sunMu)*(STM+0.5*time.simDt*k2stm);
   k4stm = orbitalState.delFdelX(xvsc+time.simDt*k3sc,constant.sunMu)*(STM+time.simDt*k3stm);
   STM = STM + time.simDt/6*(k1stm+2*k2stm+2*k3stm+k4stm);
   %% STMを用いてリファレンスの状態量と誤差共分散行列の更新
   X_bar = STM * X_hat;
   P_bar = STM * P * STM.';
   % リファレンスの状態量の時の観測量を計算する
   Y_bar = observe.calcG(X_bar,observe.xve(:,i_2),observe.xvg(:,i_2),constant);
   y     = Y - Y_bar;
   H_childa = observe.delGdelX(X_bar,observe.xve(:,i_2),observe.xvg(:,i_2),constant);
   %% カルマンゲインの計算と推定値の更新
   K = P_bar * H_childa.'/(H_childa*P_bar*H_childa.' + R);
   X_hat = X_bar + K*y;
   P = (eye(7) - K*H_childa)*P_bar;  
   % 記録する
   scEst.clockErrorCorrection(i_2)= X_hat(1);
   scEst.pos(:,i_2) = X_hat(2:4);
   scEst.vel(:,i_2) = X_hat(5:7);
end

%plot
% 時計誤差の時間履歴
figure(1)
title('clock error of spacecraft')
hold on
plot(time.list-time.list(1), error.clock - scEst.clockErrorCorrection)
xlabel('time [s]')
ylabel('clock error [s]')
hold off

% 探査機位置誤差の時間履歴
figure(2)
tiledlayout(2,1)
nexttile
title('position error of estimated value')
hold on
plot(time.list-time.list(1), scEst.pos(1,:) - scTrue.pos(1,:))
plot(time.list-time.list(1), scEst.pos(2,:) - scTrue.pos(2,:))
plot(time.list-time.list(1), scEst.pos(3,:) - scTrue.pos(3,:))
xlabel('time [s]')
ylabel('position error [m]')
legend('x', 'y', 'z')
ylim([-1000 1000])
hold off
nexttile
title('position error of initail guess')
hold on
plot(time.list-time.list(1), scInitialGuess.pos(1,:) - scTrue.pos(1,:))
plot(time.list-time.list(1), scInitialGuess.pos(2,:) - scTrue.pos(2,:))
plot(time.list-time.list(1), scInitialGuess.pos(3,:) - scTrue.pos(3,:))
xlabel('time [s]')
ylabel('position error [m]')
legend('x', 'y', 'z')
ylim([-1000 1000])
hold off

% 探査機位置の真値と推定値の比較
figure(3)
title('estemated/true position of spacecraft')
hold on
plot3(scTrue.pos(1,:),scTrue.pos(2,:),scTrue.pos(3,:))
plot3(scEst.pos(1,:),scEst.pos(2,:),scEst.pos(3,:))
plot3(scInitialGuess.pos(1,:),scInitialGuess.pos(2,:),scInitialGuess.pos(3,:))
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
legend('true', 'estimated', 'initial guess')
hold off



% hold on
% plot(time.list, ephemData.earth(1,:))
% plot(time.list - observe.ltd, observe.xve(1,:))
% hold off