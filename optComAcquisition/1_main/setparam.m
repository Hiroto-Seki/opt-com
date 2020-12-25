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
    
    
    %% parameter related to time
    % simulation timeStep[s]
    time.simDt = 10;
    % number of time step
    time.stepNum = 3000; 
    % simulateion start time (ephemeris time)
    time.t0 = cspice_str2et('2030/01/01 00:00:00 UTC');
    time.t0Ephemeris = 0;
    % 基準となる時刻のリスト
    time.list = linspace(time.t0,time.t0+time.simDt * time.stepNum,time.stepNum+1);
    
    %% parameter related to error
    % 初期時計誤差
    error.clockSigma = 1e-1; %初期時計誤差(秒), 100ppbで，約2ヶ月分蓄積した場合
    error.clock0     = error.clockSigma * randn;
    error.randomClock        = 1e-7;    %ランダム時計誤差. 帯域幅に相当
    % 初期宇宙機軌道誤差[km]. (1軸あたりの誤差は1/√3 になる)
    error.scPosSigma = 1e5; %変更した 
    % 適当に0.1km/s程度の誤差とする
    error.scVelSigma = 1e0; %変更した
    % ダイナミクスの不確定性の標準偏差(探査機)
    error.dynamics = 1e-10;
    % STTの精度
    error.stt = 10 * 10^-6 ; %ISSL unit→10urad(cross bore sight), ASTRO APS(2kg,5W程度)→1arcsec以下(4.85urad)
    % 参考: https://blog.satsearch.co/2019-11-26-star-trackers-the-cutting-edge-celestial-navigation-products-available-on-the-global-space-marketplace
    % 加速度センサの精度. 擾乱とかの方が大きいかも・・
    error.accel = 1e-12; %ちょっとサイズが大きいけど https://www.researchgate.net/publication/268554054_High-performance_Accelerometer_for_On-orbit_Spacecraft_Autonomy  
    % duration time(探査機が光を受けて返すまでの時間)の誤差
    error.duration = 1e-10;                    % 高精度にできると仮定
    % 地上局のポインティング精度
    error.gsPoint = 20*1e-9; %S-340 Piezo Tip/Tilt-Mirror Platform: High-Dynamics for Optics to 100 mm (4") Dia. mirrorが小さく，resolution=制御精度ではないが．．．20nrad
    
    %% ground station
    gs.lat  = 36.1325063*cspice_rpd();
    gs.lon =138.3607113*cspice_rpd();
    gs.alt  = 1.456;
    %% laser form gs to sc
    % 探索範囲(rad)
    gs.searchArea     = (error.scPosSigma/(SSD.AU*10)) * 3 ; %10AUくらいを想定．3sigmaをカバーする %上書きされる値
    gs.searchStep     = min(2e-6, ceil(gs.searchArea/20*1e8)*1e-8 ); %探索時の1stepあたりの間隔(rad) %上書きされる値
    gs.searchTimeStep = 2e-2 ;  %適当．SOTAの資料にこれに相当するかは分からないが20msの記述あり %上書きされる値
    % 探索1回にかかる時間
    time.obs = (2 * gs.searchArea/gs.searchStep)^2 * gs.searchTimeStep; %上書きされる値
    % 探索1回にかかる時間がシミュレーションの何stepに相当するか
    time.obsStep = ceil(time.obs/time.simDt); %上書きされる値
    % 地上局のレーザーに関するパラメータ
    gs.peakPower     = 370*10^3; 
    gs.tEff  = 0.7;
    gs.atmosphereEff = 0.73;
    gs.tAperture      = 1;
    gs.gamma         = 0.2; % telescope obscuration ratio
    gs.wavelength    = 1030 * 10^-9;
    gs.alpha         = 1.12 - 1.30 *  gs.gamma^2 + 2.12 * gs.gamma^4;
    gs.tAntGain      = (pi * gs.tAperture/ gs.wavelength)^2 * 2/gs.alpha * (exp(-gs.alpha^2) - exp(-gs.alpha^2*gs.gamma^2))^2;
    
    
    
    % 探査機
    sc.aperture = 0.22; % [m]
    sc.gamma    = 0.2;
    sc.rAntGain = (pi * sc.aperture/ gs.wavelength)^2  * (1- sc.gamma^2)^2;
    sc.rEff = 0.7;
    sc.fL   = 25 * 1e-3;   % 焦点距離
    % QD センサ(探査機)
    k = 1.381*1e-23;     %ボルツマン定数
    scT = 298;             % 絶対温度
    scRsh = 200 * 10^6;     % 並列抵抗
    sc.qdBw =  2*10^7;   % 帯域幅. 
    sc.qdFov = 200*1e-6; % 視野角 STTの姿勢決定精度が5urad(1sigma) → PAAより大きい方がいい
    sc.qdIj  = (4 * k * scT * sc.qdBw/scRsh)^0.5; %熱雑音電流
    sc.qdS   = 0.6;     % 受光感度[A/W]
    sc.qdId  = 1.5*1e-10;   % 暗電流
    sc.qdPhi = sc.fL * sc.qdFov; % 受光スポット径
    
    %% laser from sc to gs
    % European Deep Space Optical Communication Programを参考にした
    % 探査機
    sc.meanPower = 10;
    sc.peakPower = 640; 
    sc.tEff = 0.7;
    sc.atmosphereEff = 0.73;
    sc.wavelength = 1550 * 10^-9;
    sc.alpha = 1.12 - 1.30 *  sc.gamma^2 + 2.12 * sc.gamma^4;
    sc.tAntGain      = (pi * sc.aperture/ sc.wavelength)^2 * 2/sc.alpha * (exp(-sc.alpha^2) - exp(-sc.alpha^2*sc.gamma^2))^2;
    
    % 地上局
    gs.rAperture = 10; 
    gs.rAntGain = (pi * gs.rAperture/ sc.wavelength)^2  * (1- gs.gamma^2)^2;
    gs.rEff = 0.5;
    % QD センサに相当するもの(地上局)
    gsT = 1;           % 絶対温度
    gsRsh = 50 * 10^6;     % 並列抵抗
    gs.qdBw =  2.5e8;   % 帯域幅. pulse widthを大きくできれば可能 2e7に設定．仮置きで小さな値を設定
    gs.qdFov = 200*1e-6  ;  % 視野角 STTの姿勢決定精度が5urad(1sigma) 
    gs.qdIj  = (4 * k * gsT * sc.qdBw/gsRsh)^0.5; %熱雑音電流
    gs.qdS   = 1.16;     % 受光感度[A/W]
    gs.qdId  = 1.6*1e-16;   % 暗電流 dark count=1[kcps]に相当．Large-Area 64-pixel Array of WSi Superconducting　Nanowire Single Photon Detectors　を参考
        
    %% 地上局, 宇宙機の位置の設定  
    % 地上局
    gsTrue = GroundStation(gs,constant,time);
    % 地球
    earth  = CelestialBody(time,"Earth");
    earth.getEphem(time);
    % 宇宙機
    sc.state0 = cspice_spkezr('699', time.t0Ephemeris + time.t0,'ECLIPJ2000', 'NONE', '10'); %推定値の初期値. とりあえず土星にしている
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
    scEstByScSeq.R1wSc = [error.stt^2*eye(2),                                     zeros(2,6);         % 測角 (受信電力で書き換える)
                              zeros(3,2),    1e0*        error.accel^2*eye(3),zeros(3,3);         % 加速度計
                              zeros(2,5),       (gs.searchStep^2+error.gsPoint^2)*eye(2),zeros(2,1);         % uplinkの送信方向
                              zeros(1,7),            (1e0*error.randomClock * constant.lightSpeed)^2];          %1wayの測距 おそらくどこかで桁落ち誤差が発生しているので精度を落としている
    scEstByScSeq.R2wSc = [scEstByScSeq.R1wSc, zeros(8,1); zeros(1,8),  max((1e0 * error.randomClock * constant.lightSpeed)^2, 1e3^2)];% 2wayの測距  
                          
    scEstByGsSeq.X             = [0;sc.state0];
    scEstByGsSeq.P             = [error.clockSigma^2,                                               zeros(1,6);
                                          zeros(3,1), 1/3 * error.scPosSigma^2 * eye(3),                  zeros(3,3);
                                                                       zeros(3,4), 1/3 * error.scVelSigma^2 * eye(3)];
    scEstByGsSeq.P_list        = zeros(7,7,length(time.list));
    scEstByGsSeq.P_list(:,:,1) = scEstByGsSeq.P;
    scEstByGsSeq.R2wGs = [error.stt^2*eye(2),                                           zeros(2,5);         % 測角 (受信電力で書き換える)
                              zeros(3,2), 1e0* error.accel^2*eye(3),                          zeros(3,2);         % 加速度計 
                              zeros(1,5),     (1e0*error.randomClock * constant.lightSpeed)^2,                                    0;
                              zeros(1,6),     (1e0*error.randomClock * constant.lightSpeed)^2];         % 2wayの測距  

                          
    % EKFに使用するパラメーター
    ekf.sigmaN = 2;
    
    
    % UKFに使用するパラメーター
    ukf.n = 7; % 推定する状態量の数(クロックオフセット，位置3・速度3)
    ukf.alpha = 0.5; % 0~1の間
    % 棄却シグマ値レベル．どのくらいが妥当なのか？？→とりあえず全て3シグマにしておく
    ukf.sigmaN = 2;
    
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
   scTrue.ur_counter=0;
   scTrue.dt_counter=0;
   scTrue.ur2w_counter = 0;
   % 初期化
   time.lastSearch = 0;
   time.sc2wayget  = time.list(length(time.list)) + time.simDt * 2;  %大きい値で初期化しておく
   
   
end