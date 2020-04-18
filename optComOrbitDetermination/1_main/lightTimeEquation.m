% solve light time eqaution estimating the orbit of spacecraft, ground
% station

clear all; close all; clc

%% setting parameter
rng('default');
rng(1)
% 定数
constant.sunMu         = 1.32712438e11;   % km^3/s^2 太陽の重力定数
constant.earthMu       = 3.98600441e5;    % km^3/s^2 地球の重力定数
constant.saturnMu      = 3.7931187e7;     %[km/s^2]
constant.au            = 149597870.7;     % km 太陽-地球間距離
constant.lightSpeed    = 299792.458;      % [km/s] 光速 
constant.earthRadius   = 6.371e3;         % [km] 地球半径
constant.earthRotation = 2*pi/24/60/60;   % [rad/s] 地球自転速度
constant.eathAxis      = 23.4/180*pi;     % 地軸の傾き 

%% 時間に関するパラメーター
% simulation timeStep[s]
time.simDt = 100;
% 基準時刻は，time.listの中で何番目か
time.refId = 600; 
% 基準時刻0としてどの範囲の時刻を計算するか(-6000秒~6000秒をデフォルトとする)
time.list = linspace((-time.refId+1)*time.simDt, (time.refId)*time.simDt,time.refId*2);

% 基準時刻の設定
time.t0 = juliandate(2030,1,1,0,0,0);

%% 誤差に関するパラメーター
% 初期時計誤差
error.clock = 10   * randn; %時計誤差(秒)
error.time0 = error.clock/ (60*60*24); %時計誤差(日)
% 初期探査機軌道誤差[km]. 地球からみて1μradくらいの誤差にする
% error.scPos0 = randn(3,1) *  150;
error.scPos0 = randn(3,1) *  0;
% 適当に0.1km/s程度の誤差とする
% error.scVel0 = randn(3,1) *  0.1;
error.scVel0 = randn(3,1) *  0;
%%  初期推定値の設定
% 地球の推定時刻のephemeris
earthEst = orbitalState('Earth','estimate','SolarSystem');
[earthEst.pos0,earthEst.vel0]=planetEphemeris(time.t0,'SolarSystem','Earth');
earthEst.pos0 =  reshape(earthEst.pos0,[3,1]);
earthEst.vel0 =  reshape(earthEst.vel0,[3,1]);
% とりあえず，探査機軌道の推定値 = 土星軌道の推定値とする
scEst = orbitalState('spacecraft','estimate','SolarSystem');
[scEst.pos0,scEst.vel0]=planetEphemeris(time.t0,'SolarSystem','Saturn');
scEst.pos0 =  reshape(scEst.pos0,[3,1]);
scEst.vel0 =  reshape(scEst.vel0,[3,1]);

% 地上局
gsEst = groundState('estimate');
% 地上局の設定(ECI座標系で)(適当)
gsEst.pos0 =  [constant.earthRadius;0;0];

%% 真値の設定
%基準時刻のephemerisの真値
earthTrue = orbitalState('Earth','true','SolarSystem');
[earthTrue.pos0,earthTrue.vel0]=planetEphemeris(time.t0,'SolarSystem','Earth');
earthTrue.pos0 =  reshape(earthTrue.pos0,[3,1]);
earthTrue.vel0 =  reshape(earthTrue.vel0,[3,1]);
% [saturnPosTrue0,saturnVelTrue0]=planetEphemeris(time.t0+error.time0,'SolarSystem','Saturn');
% 探査機の土星中心での位置速度
scTrue = orbitalState('spacecraft','true','SolarSystem');
scTrue.pos0 = scEst.pos0 + error.scPos0;
scTrue.vel0 = scEst.vel0 + error.scVel0;

% 地上局の真値
gsTrue = groundState('true');
[gsTrue.pos0,~] = groundState.earthRotation(gsEst.pos0, 0, constant);

%% 探査機は軌道伝搬によってtime.list分の位置，速度を得る
scEst.getOrbitTwoBody(time, constant);
scTrue.getOrbitTwoBody(time,constant);

%% 地球の位置，速度を伝搬して求める
earthEst.getOrbitTwoBody(time, constant);
earthTrue.getOrbitTwoBody(time,constant);

%% 地上局の位置，速度を自転の運動(t)の関数によって求める
gsEst.getTrajectoryEarthRotation(time,constant)
gsTrue.getTrajectoryEarthRotation(time,constant);

%% 軌道を時刻の関数にする(シミュレーション区間を通しての関数とする)
% 二字曲線近似 or r(t+dt)=r(t)+v(t)dt+0.5a(t)*dt^2， Chebyshev coefficients　→あとまわし

%% 観測する(観測は真値を伝搬して求めていく)
observe = observeState(time);
% 初期予想時計誤差
clockErrorCorrection = 0;
for i = 1:length(time.list)
    % まずは観測量を得る
    observe.fromSc2Gs(time,earthTrue,gsTrue,scTrue,constant,error.clock,i);
    % 観測量から時計誤差を推測する
    observe.calcClockOffset(time,earthEst,gsEst,scEst,constant,i,clockErrorCorrection)
    clockErrorCorrection = observe.clockError(i);
end

hold on
plot(time.list/60/60/24, observe.clockError)
plot(time.list/60/60/24, error.clock * ones(1, length(time.list)))
xlabel('year')
ylabel('light time delay [sec]')
legend('estimated','true')
hold off


% 地球位置の

%% 観測から伝搬遅延を求める(ここでは軌道は時間の関数としたものを用いる)→時計誤差を見積もりたい
% observe.calcClockOffset(time,earthEst,gsEst,scEst,constant)


% ここまでを今日中にやる

%% 見積もった伝搬遅延から，指向させる(初期軌道がずれていればちゃんと指向できないが，真値と推定値両方を計算して，差分を見たい)


%%% ここまではRGまでにやりたい
%% 将来的には，
% 1. 観測を増やして軌道を推定する

% 2. 推定した軌道をもとに地上局に指向する
% ここで，どのように探索するか検討(ただ平均値の方向に向くか，探すか)

% 3. 地上局で受信できる
% 観測を重ねて地上局が，探査機の軌道を推定する．(探査機の時計のずれ，軌道)
% ここでは，完全にtwo wayが確立できたと仮定し，探査機の時計のずれは完全にわかったとする

% 4. 地上局のアップリンクによって，探査機は時計(地上局と一致)と軌道を補正する．

% 5 最終的な指向精度を見積もる．



