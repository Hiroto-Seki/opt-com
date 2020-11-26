function [constant,time,error,gs,sc,gsTrue,earth,scTrue,scEstByScEkf,scEstByGsEkf] = setparam(SSD)
    %% constant value
    constant.sunMu         = SSD.GM(10);      % km^3/s^2 gravity constant of the sun
    constant.earthMu       = SSD.GM(399);     % km^3/s^2 gravity constant of the 
    constant.saturnMu      = SSD.GM(799);     %[km/s^2]
    constant.lightSpeed    = 299792.458 ;      % [km/s] light speed 
    constant.earthRadius   = SSD.r(399);      % [km] earth radius
    constant.earthRotation = 2*pi/24/60/60;   % [rad/s] earth rotation speed
    constant.eathAxis      = 23.4/180*pi;     % inclination of the earth axis
    constant.earthF = 1/298.257223560; % https://topex.ucsd.edu/geodynamics/14gravity1_2.pdf
    constant.elementaryCharge = 1.602 * 1e-19;
    
    %% parameter related to time
    % simulation timeStep[s]
    time.simDt = 10;
    % number of time step
    time.stepNum = 2160; 
    % simulateion start time (ephemeris time)
    time.t0 = cspice_str2et('2030/01/01 00:00:00 UTC');
    % 基準となる時刻のリスト
    time.list = linspace(time.t0,time.t0+time.simDt * time.stepNum,time.stepNum+1);
    
    %% parameter related to error
    % 初期時計誤差
    error.clockSigma = 0.5; %初期時計誤差(秒), 100ppbで，約2ヶ月分蓄積した場合
    error.clock0     = error.clockSigma * randn;
    error.randomClock        = 1e-8;         %ランダム時計誤差
    % 初期宇宙機軌道誤差[km]. (1軸あたりの誤差は1/√3 になる)
    error.scPosSigma = 2e4;
    % 適当に0.1km/s程度の誤差とする
    error.scVelSigma = 5e-1;
    % ダイナミクスの不確定性の標準偏差(探査機)
    error.dynamics = 1e-10;
    % STTの精度
    error.stt = 3 * 10^-6; %1urad
    % 加速度センサの精度
    error.accel = 1e-10;   
    % duration time(探査機が光を受けて返すまでの時間)の誤差
    error.duration = 1e-10;                    % 高精度にできると仮定
    
    %% ground station
    gs.lat  = 36.1325063*cspice_rpd();
    gs.lon =138.3607113*cspice_rpd();
    gs.alt  = 1.456;
    %% laser form gs to sc
    % 探索範囲(rad)
    gs.searchArea     = (error.scPosSigma/(SSD.AU*10)) * 3 ; %10AUくらいを想定．3sigmaをカバーする
%     gs.searchArea     = 2e-5;
    gs.searchStep     = 2e-6; %探索時の1stepあたりの間隔(rad)
    gs.searchTimeStep = 2e-2; %適当．SOTAの資料にこれに相当するかは分からないが20msの記述あり
    % 探索1回にかかる時間
    time.obs = (2 * gs.searchArea/gs.searchStep)^2 * gs.searchTimeStep;
    % 探索1回にかかる時間がシミュレーションの何stepに相当するか
    time.obsStep = ceil(time.obs/time.simDt);
    % 地上局のレーザーに関するパラメータ
    gs.peakPower     = 370*10^3 ;
%   gs.meanPower     = 2.4*10^3;
    gs.tEff  = 0.7;
    gs.atmosphereEff = 0.73;
    gs.tAperture      = 1;
    gs.gamma         = 0.2; % telescope obscuration ratio
    gs.wavelength    = 1030 * 10^-9;
    gs.alpha         = 1.12 - 1.30 *  gs.gamma^2 + 2.12 * gs.gamma^4;
    gs.tAntGain      = (pi * gs.tAperture/ gs.wavelength)^2 * 2/gs.alpha * (exp(-gs.alpha^2) - exp(-gs.alpha^2*gs.gamma^2))^2;
    
    % 探査機
    sc.aperture = 0.22; 
    sc.gamma    = 0.2;
    sc.rAntGain = (pi * sc.aperture/ gs.wavelength)^2  * (1- sc.gamma^2)^2;
    sc.rEff = 0.7;
    % QD センサ(探査機)
    k = 1.381*1e-23;     %ボルツマン定数
    scT = 298;             % 絶対温度
    scRsh = 200 * 10^6;     % 並列抵抗
    sc.qdBw =  2*10^7;   % 帯域幅. 
    sc.qdFov = 200*1e-6; % 視野角 STTの姿勢決定精度が5urad(1sigma) → PAAより大きい方がいいかも 250 uradにしたい...
    sc.qdIj  = (4 * k * scT * sc.qdBw/scRsh)^0.5; %熱雑音電流
    sc.qdS   = 0.6;     % 受光感度[A/W]
    sc.qdId  = 1.5*1e-10;   % 暗電流
    
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
    gs.qdFov = 200*1e-6;  % 視野角 STTの姿勢決定精度が5urad(1sigma) 
    gs.qdIj  = (4 * k * gsT * sc.qdBw/gsRsh)^0.5; %熱雑音電流
    gs.qdS   = 1.16;     % 受光感度[A/W]
    gs.qdId  = 1.6*1e-16;   % 暗電流 dark count=1[kcps]に相当．Large-Area 64-pixel Array of WSi Superconducting　Nanowire Single Photon Detectors　を参考
        
    %% 地上局, 宇宙機の位置の設定  
    % 地上局
    gsTrue = GroundStation(gs,constant,time);
    % 地球
    earth  = CelestialBody(time,constant.sunMu,"Earth");
    earth.getEphem();
    % 宇宙機
    sc.state0 = cspice_spkezr('699', time.t0,'ECLIPJ2000', 'NONE', '10'); %推定値の初期値. とりあえず土星にしている
    % 真値(誤差を入れる)
    scTrue                  = Spacecraft(time,constant.sunMu);
%     scTrue.state            = sc.state0 + [ 1/3^0.5 * randn(3,1) * error.scPosSigma ; 1/3^0.5 * randn(3,1) * error.scVelSigma ];
      scTrue.state            = sc.state0    +[23623.3208393106;-15169.0659456738;-22192.2607700304;-0.422775620003898;-0.286332433228975;-0.279340382236986 ];
    % 宇宙機自身がEKFで推定した値
    scEstByScEkf            = Spacecraft(time,constant.sunMu);
    scEstByScEkf.state      = sc.state0;
    scEstByScEkf.clockError = error.clock0; 
    % 地上局がEKFで推定した値
    scEstByGsEkf            = Spacecraft(time,constant.sunMu);
    scEstByGsEkf.state      = sc.state0;
    scEstByGsEkf.clockError = error.clock0;     
    
    %% 宇宙機の姿勢を表す. (推定値)
    % 慣性空間上での宇宙機から見た地上局方向. 
    sc2gs = earth.state(:,1) + gsTrue.state(:,1) - scEstByScEkf.state(:,1);
    % 姿勢は適当．だいたい地上局方向に光学系の座標系が向くように設定する.(厳密には一致させなくても良い)
    roll  = pi/2;
    pitch = atan2(sc2gs(2),sc2gs(1));
    yaw   = 0;
    % 宇宙機の姿勢の初期値
    scTrue.attState = [roll + randn * error.stt ; pitch + randn * error.stt ;yaw + randn * error.stt];

    %% 推定値 (clockのオフセット + 宇宙機の位置・速度)
    scEstByScEkf.X             = [0;sc.state0];
    scEstByScEkf.P             = zeros(7);
    scEstByScEkf.P(1,1)        = error.clockSigma^2;
    scEstByScEkf.P(2,2)        = error.scPosSigma^2;
    scEstByScEkf.P(3,3)        = error.scPosSigma^2;
    scEstByScEkf.P(4,4)        = error.scPosSigma^2;
    scEstByScEkf.P(5,5)        = error.scVelSigma^2;
    scEstByScEkf.P(6,6)        = error.scVelSigma^2;
    scEstByScEkf.P(7,7)        = error.scVelSigma^2;
    scEstByScEkf.P_list        = zeros(7,7,length(time.list));
    scEstByScEkf.P_list(:,:,1) = scEstByScEkf.P;
    scEstByScEkf.R = [error.stt^2,0,0,0,0,0;         % 測角 (受信電力で書き換える)
                      0,error.stt^2,0,0,0,0;         % 測角 (受信電力で書き換える)
                      0,0,error.accel^2,0,0,0;       % 加速度計
                      0,0,0,error.accel^2,0,0;       % 加速度計
                      0,0,0,0,error.accel^2,0;       % 加速度計
                      0,0,0,0,0,error.randomClock^2];  %1wayの測距  %2wayも得られたら，7x7の行列に拡張する    

    scEstByGsEkf.X             = [0;sc.state0];
    scEstByGsEkf.P             = zeros(7);
    scEstByGsEkf.P(1,1)        = error.clockSigma^2;
    scEstByGsEkf.P(2,2)        = error.scPosSigma^2;
    scEstByGsEkf.P(3,3)        = error.scPosSigma^2;
    scEstByGsEkf.P(4,4)        = error.scPosSigma^2;
    scEstByGsEkf.P(5,5)        = error.scVelSigma^2;
    scEstByGsEkf.P(6,6)        = error.scVelSigma^2;
    scEstByGsEkf.P(7,7)        = error.scVelSigma^2;
    scEstByGsEkf.P_list        = zeros(7,7,length(time.list));
    scEstByGsEkf.P_list(:,:,1) = scEstByScEkf.P;
    scEstByGsEkf.R = [error.stt^2,0,0,0,0,0,0;         % 測角 (受信電力で書き換える)
                      0,error.stt^2,0,0,0,0,0;         % 測角 (受信電力で書き換える)
                      0,0,error.accel^2,0,0,0,0;       % 加速度計
                      0,0,0,error.accel^2,0,0,0;       % 加速度計
                      0,0,0,0,error.accel^2,0,0;       % 加速度計
                      0,0,0,0,0,error.randomClock^2,0; %1wayの測距  
                      0,0,0,0,0,0,error.randomClock^2];%2wayの測距  

                  
   % 送受信した回数の初期化
   gsTrue.ut_counter=0;
   gsTrue.dr_counter=0;
   scTrue.ur_counter=0;
   scTrue.dt_counter=0;
   % 初期化
   time.lastSearch = 0;
   
end