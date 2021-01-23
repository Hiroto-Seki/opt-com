function [constant,time,error,gs,sc,gsTrue,earth,scTrue,scEstByScSeq,scEstByGsSeq,ekf,ukf] = setparam(SSD)
    %% constant value
    constant.sunMu         = SSD.GM(10);      % km^3/s^2 gravity constant of the sun
    constant.earthMu       = SSD.GM(399);     % km^3/s^2 gravity constant of the 
    constant.saturnMu      = SSD.GM(799);     %[km/s^2]
    constant.lightSpeed    = 299792.458;     % [km/s] light speed 
    constant.earthRadius   = SSD.r(399);      % [km] earth radius
    constant.earthRotation = 2*pi/24/60/60;   % [rad/s] earth rotation speed
    constant.eathAxis      = 23.4/180*pi;     % inclination of the earth axis
    constant.earthF = 1/298.257223560; % https://topex.ucsd.edu/geodynamics/14gravity1_2.pdf
    constant.elementaryCharge = 1.602 * 1e-19;
    constant.boltzmann = 1.381*1e-23;         %ボルツマン定数
    
    
    %% parameter related to time
    % simulation timeStep[s]
    time.simDt = 60;
    % number of time step
    time.stepNum = 10000; 
    % simulateion start time (ephemeris time)
%     time.t0 = cspice_str2et('2030/01/01 00:00:00 UTC');
    time.t0 = cspice_str2et('2034 DEC 30 20:55:10.305362');
    time.t0Ephemeris = 0;
    % LOSの何日前から軌道決定するか(とりあえず1週間前から計算)
    time2LOS = 24 * 60 * 60 * 7;
    time.t0  = time.t0 - time2LOS;
    
    % 基準となる時刻のリスト
    time.list = linspace(time.t0,time.t0+time.simDt * time.stepNum,time.stepNum+1);
    
    % 宇宙機の初期値
    % sc.state0 = cspice_spkezr('699', time.t0Ephemeris + time.t0,'ECLIPJ2000', 'NONE', '10'); %推定値の初期値. とりあえず土星にしている
    sc.state0 = [-738505094.095115; 1108698819.51426; 15674510.6716146; -5.19045777910223;-0.349637820348331; -0.212238899522963];
    % 宇宙機の初期値をtime2LOSだけ逆伝搬して求める
    sc.state0 = Spacecraft.timeUpdate_sc(sc.state0,constant.sunMu,-time2LOS, time.simDt);
    
    
    %% parameter related to error
    % 初期時計誤差
    error.clockSigma = 1e-1; %初期時計誤差(秒), 100ppbで，約2ヶ月分蓄積した場合
    error.clock0     = error.clockSigma * randn;
    error.randomClock        = 1e-7;    %ランダム時計誤差. 帯域幅に相当
    % 初期宇宙機軌道誤差[km]. (1軸あたりの誤差は1/√3 になる)
    error.scPosSigma = 1e5; %変更した 
    % 適当に0.1km/s程度の誤差とする
    error.scVelSigma = 5e-1; %変更した
    % ダイナミクスの不確定性の標準偏差(探査機) 太陽輻射厚が100kg,
    % 10m^2で，4.6e-12km/s^2で，それより少し小さめの値に設定した
    error.dynamics = 1e-11;
    % STTの精度
    error.stt = 4.85 * 10^-6 * 0.364 ; %ISSL unit→10urad(cross bore sight), ASTRO APS(2kg,5W程度)→1arcsec以下(4.85urad)
    % 参考: https://blog.satsearch.co/2019-11-26-star-trackers-the-cutting-edge-celestial-navigation-products-available-on-the-global-space-marketplace
    % 加速度センサの精度. 擾乱とかの方が大きいかも・・
    error.accel = 1e-12; %ちょっとサイズが大きいけど https://www.researchgate.net/publication/268554054_High-performance_Accelerometer_for_On-orbit_Spacecraft_Autonomy  
%     % duration time(探査機が光を受けて返すまでの時間)の誤差
%     error.duration = 1e-10;                    % 高精度にできると仮定
    % 地上局のポインティング精度. 
    error.gsPoint = 1*1e-7; %S-340 Piezo Tip/Tilt-Mirror Platform: High-Dynamics for Optics to 100 mm (4") Dia. mirrorが小さく，resolution=制御精度ではないが．．．20nrad
    % 宇宙機のポインティング精度
    error.scPoint = 1*1e-7;
    
    %% ground station
    gs.lat  = 36.1325063*cspice_rpd();
    gs.lon = 138.3607113*cspice_rpd();
    gs.alt  = 1.456;
    %% laser form gs to sc
    % 探索範囲(rad)
    gs.searchArea     = (error.scPosSigma/(SSD.AU*10)) * 3 ; %10AUくらいを想定．3sigmaをカバーする %上書きされる値
    gs.searchStep     = min(0.8e-6, ceil(gs.searchArea/20*1e8)*1e-8 ); %探索時の1stepあたりの間隔(rad) %上書きされる値
    gs.switchTime     = 2e-2; %適当．SOTAの資料にこれに相当するかは分からないが20msの記述あり
    gs.searchTimeStep = gs.switchTime + 0 ;   % 送信する情報量によって上書きする
    % 探索1回にかかる時間
    time.obs = (2 * gs.searchArea/gs.searchStep)^2 * gs.searchTimeStep; %上書きされる値
    % 探索1回にかかる時間がシミュレーションの何stepに相当するか
    time.obsStep = ceil(time.obs/time.simDt); %上書きされる値
    % 地上局のレーザーに関するパラメータ
    gs.peakPower     = 370*10^3;
    gs.peakWidth     = 128 * 1e-9;
    gs.Ppm_up        = 256;
    gs.symbolRate_up = 0.66;
    gs.tEff          = 0.7 * 0.7; %バンドパス損失も含めた
    gs.atmosphereEff = 0.73;
    gs.tAperture     = 1;
    gs.gamma         = 0.2; % telescope obscuration ratio
    gs.wavelength_up = 1030 * 10^-9;
    gs.alpha         = 1.12 - 1.30 *  gs.gamma^2 + 2.12 * gs.gamma^4;
    gs.tAntGain      = (pi * gs.tAperture/ gs.wavelength_up)^2 * 2/gs.alpha * (exp(-gs.alpha^2) - exp(-gs.alpha^2*gs.gamma^2))^2; 
    
    % 探査機
    sc.aperture = 0.22; % [m]
    sc.gamma    = 0.2;
    sc.rAntGain = (pi * sc.aperture/ gs.wavelength_up)^2  * (1- sc.gamma^2)^2;
    sc.rEff = 0.7 * 0.7; % バンドパス損失も含めた
    sc.fL   = 25 * 1e-3;   % 焦点距離
    % QD センサ(探査機)
    k = constant.boltzmann;     %ボルツマン定数
    sc.T = 263;             % 絶対温度
    sc.Rsh = 200 * 10^6 * 10^(-0.0375*(sc.T-298));     % 並列抵抗
    sc.qdBw =  16.08*1e6;   % 帯域幅. 
    sc.qdGain = 1.25;      %フォトダイオードのゲイン
    sc.qdF    = 2 + sc.qdGain * 0.028; %過剰ノイズ比
    sc.qdFov = 385*1e-6; % 視野角 STTの姿勢決定精度が5urad(1sigma) → PAAより大きい方がいい
    sc.qdIj  = (4 * k * sc.T * sc.qdBw/sc.Rsh)^0.5; %熱雑音電流
    sc.qdS   = 0.95;     % 受光感度[A/W]
    sc.qdId  = 0.15 * 1e-9 * 1.09^(sc.T-298) ;   % 暗電流
    sc.qdPhi = sc.fL * sc.qdFov; % 受光スポット径
    sc.reqBer_up = 1e-5; %要求ビット誤り率
    
    % 計算されるパラメーター
    % PPM変調の1slotあたりの平均電力
    if 2/sc.qdBw > gs.peakWidth
        gs.slotPower = gs.peakPower * (gs.peakWidth/(2/sc.qdBw));
    else
        gs.slotPower = gs.peakPower;
    end
    % uplinkに要求されるS/N比(slotで評価した値)
    sc.reqSnr_up = 16 * erfcinv(2*sc.reqBer_up)^2/(gs.symbolRate_up*log2(gs.Ppm_up));
    % 通信速度
    gs.bps_up    =  gs.symbolRate_up * log2(gs.Ppm_up)/gs.Ppm_up/(2/sc.qdBw);
    
    %% laser from sc to gs
    % European Deep Space Optical Communication Programを参考にした
    % 探査機
    sc.peakPower = 640 * 9.58;
    sc.peakWidth = 8 * 1e-9 * 9.58;
    sc.Ppm_down         = 4096;
    sc.symbolRate_down  = 0.66;
    sc.tEff = 0.7;
    sc.atmosphereEff = 0.73;
    sc.wavelength_down = 1550 * 10^-9;
    sc.alpha = 1.12 - 1.30 *  sc.gamma^2 + 2.12 * sc.gamma^4;
    sc.tAntGain      = (pi * sc.aperture/ sc.wavelength_down)^2 * 2/sc.alpha * (exp(-sc.alpha^2) - exp(-sc.alpha^2*sc.gamma^2))^2;
    
    % 地上局
    gs.rAperture = 12; 
    gs.rAntGain = (pi * gs.rAperture/ sc.wavelength_down)^2  * (1- gs.gamma^2)^2;
    gs.rEff = 0.5;
    % QD センサに相当するもの(地上局)
    gs.T = 233;           % 絶対温度
    gs.Rsh = 200 * 10^6 * 10^(-0.0375*(gs.T-298));     % 並列抵抗
    gs.qdBw =  18.68 * 1e6;   % 帯域幅. 
    gs.qdGain = 1;      %フォトダイオードのゲイン
    gs.qdF    = 2 + gs.qdGain * 0.028; %過剰ノイズ比
    gs.qdFov = 385*1e-6  ;  % 視野角 STTの姿勢決定精度が5urad(1sigma) 
    gs.qdIj  = (4 * k * gs.T * sc.qdBw/gs.Rsh)^0.5; %熱雑音電流
    gs.qdS = 0.95;
    gs.qdId  = 0.15 * 1e-9 * 1.09^(gs.T-298) ;   % Ta = 25度，Id(T) = Id(Ta) * beta^(T-ta)
    gs.reqBer_down = 1e-5; %要求ビット誤り率
    
    % 計算されるパラメーター
    % PPM変調の1slotあたりの平均電力
    if 2/gs.qdBw > sc.peakWidth
        sc.slotPower = sc.peakPower * (sc.peakWidth/(2/gs.qdBw));
    else
        sc.slotPower = sc.peakPower;
    end
    % downlinkに要求されるS/N比(slotで評価した値)
    gs.reqSnr_down = 16 * erfcinv(2*gs.reqBer_down)^2/(sc.symbolRate_down*log2(sc.Ppm_down));
    % 通信速度
    sc.bps_down    =  sc.symbolRate_down * log2(sc.Ppm_down)/sc.Ppm_down/(2/gs.qdBw); 
    
    %% 地上局, 宇宙機の位置の設定  
    % 地上局
    gsTrue = GroundStation(gs,constant,time);
    % 地球
    earth  = CelestialBody(time,"Earth");
    earth.getEphem(time);

    
    % 真値(誤差を入れる)
    scTrue                  = Spacecraft(time);
    scTrue.state            = sc.state0 + [ 1/3^0.5 * randn(3,1) * error.scPosSigma ; 1/3^0.5 * randn(3,1) * error.scVelSigma ];
    % 宇宙機自身がEKFで推定した値
    scEstByScSeq            = Spacecraft(time);
    scEstByScSeq.state      = sc.state0;
    scEstByScSeq.clockError = error.clock0; 
    % 地上局がEKFで推定した値
    scEstByGsSeq            = Spacecraft(time);
    scEstByGsSeq.state      = sc.state0;
    scEstByGsSeq.clockError = error.clock0;     

    %% 推定値 (clockのオフセット + 宇宙機の位置・速度)
    scEstByScSeq.X             = [0;sc.state0];
    scEstByScSeq.P             = [error.clockSigma^2,                                               zeros(1,6);
                                          zeros(3,1), 1/3 * error.scPosSigma^2 * eye(3),                  zeros(3,3);
                                                                       zeros(3,4), 1/3* error.scVelSigma^2 * eye(3)];
    scEstByScSeq.P_list        = zeros(7,7,length(time.list));
    scEstByScSeq.P_list(:,:,1) = scEstByScSeq.P;
    
    
   % 使う観測の設定0or1で記述する. 0は使用しない. 1は使用する
   scEstByScSeq.useObs.direction_ur =1;  %uplinkを宇宙機が受信する角度
   scEstByScSeq.useObs.direction_ut =1;  %uplinkを地上局が送信する角度 　0にする.1にすると，たぶんうまくいけば軌道決定精度がかなり上がるが，うまくデバッグできなかった
   scEstByScSeq.useObs.accel_ur =0;      %uplinkを宇宙機が受信する時の加速度
   scEstByScSeq.useObs.length1w_ur =1;   %地上局→宇宙機の1way測距
   scEstByScSeq.useObs.length2w_ur =1;   %宇宙機→地上局→宇宙機の2way測距
   scEstByScSeq.useObs.direction_dr =0;  %downlinkを地上局が受信する角度 
   scEstByScSeq.useObs.length1w_dr =0;   %宇宙機→地上局の1way測距       (0になる)
   scEstByScSeq.useObs.length2w_dr =0;   %地上局→宇宙機→地上局の1way測距 (0になる)
   
  % 使う観測の設定0or1で記述する. 0は使用しない. 1は使用する
   scEstByGsSeq.useObs.direction_ur =1;  %uplinkを宇宙機が受信する角度
   scEstByGsSeq.useObs.direction_ut =1;  %uplinkを地上局が送信する角度 0にする.1にすると，たぶんうまくいけば軌道決定精度がかなり上がるが，うまくデバッグできなかった
   scEstByGsSeq.useObs.accel_ur =0;      %uplinkを宇宙機が受信する時の加速度
   scEstByGsSeq.useObs.length1w_ur =0;   %地上局→宇宙機の1way測距        (0になる)
   scEstByGsSeq.useObs.length2w_ur =0;   %宇宙機→地上局→宇宙機の2way測距  (0になる)
   scEstByGsSeq.useObs.direction_dr =0;  %downlinkを地上局が受信する角度
   scEstByGsSeq.useObs.length1w_dr =1;   %宇宙機→地上局の1way測距
   scEstByGsSeq.useObs.length2w_dr =1;   %地上局→宇宙機→地上局の2way測距
    
   scEstByScSeq.R.direction_ur = error.stt^2;
   scEstByScSeq.R.direction_ut = (gs.searchStep^2);
   scEstByScSeq.R.accel_ur     = error.accel^2;
   scEstByScSeq.R.length1w_ur  = (1e0 * error.randomClock * constant.lightSpeed)^2;
   scEstByScSeq.R.length2w_ur  = (1e0 * error.randomClock * constant.lightSpeed)^2;
   scEstByScSeq.R.direction_dr = error.stt^2;  %downlinkを地上局が受信する角度
   scEstByScSeq.R.length1w_dr  = (1e0 * error.randomClock * constant.lightSpeed)^2;   %宇宙機→地上局の1way測距
   scEstByScSeq.R.length2w_dr  = (1e0 * error.randomClock * constant.lightSpeed)^2;   %地上局→宇宙機→地上局の1way測距  
                             
   scEstByGsSeq.X             = scEstByScSeq.X;
   scEstByGsSeq.P             = scEstByScSeq.P;
   scEstByGsSeq.P_list        = zeros(7,7,length(time.list));
   scEstByGsSeq.P_list(:,:,1) = scEstByGsSeq.P;
   scEstByGsSeq.R             = scEstByScSeq.R;

                          
    % EKFに使用するパラメーター
    ekf.sigmaN = 3;
    
    % UKFに使用するパラメーター
    ukf.n = 7; % 推定する状態量の数(クロックオフセット，位置3・速度3)
    ukf.alpha = 0.5; % 0~1の間
    % 棄却シグマ値レベル．どのくらいが妥当なのか？？→とりあえず全て3シグマにしておく
    ukf.sigmaN = 3;
    
    ukf.kappa  = 3 - ukf.n;
    ukf.lambda = ukf.alpha^2 * (ukf.n + ukf.kappa) - ukf.n;
    ukf.beta   = 2;
    ukf.w0_m   = ukf.lambda/(ukf.n + ukf.lambda);
    ukf.wi_m   = 1/(2*(ukf.n + ukf.lambda));
    ukf.w0_c   = ukf.lambda/(ukf.n + ukf.lambda) + (1- ukf.alpha^2 + ukf.beta);
    ukf.wi_c   = 1/(2*(ukf.n + ukf.lambda));  
                                                  
   % 送受信した回数の初期化
   gsTrue.ut_counter=0;
   gsTrue.dr_counter=0;
   gsTrue.ut2w_counter = 0;
   scTrue.ur_counter=0;
   scTrue.dt_counter=0;
   scTrue.ur2w_counter = 0;
   % 初期化
   time.lastSearch = 0;
   % 結果を格納する
   gsTrue.t_drList = [];
   gsTrue.snr_drList = [];
   scTrue.targetError_dtList = [];
   scTrue.pointError_dtList = [];
end