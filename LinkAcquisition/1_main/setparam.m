function [constant,time,error,gs,sc,X_hat,P,P_list,R] = setparam(SSD)
    %% constant value
    constant.sunMu         = SSD.GM(10);      % km^3/s^2 gravity constant of the sun
    constant.earthMu       = SSD.GM(399);     % km^3/s^2 gravity constant of the 
    constant.saturnMu      = SSD.GM(799);     %[km/s^2]
    constant.lightSpeed    = 299792.458;      % [km/s] light speed 
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
    time.list = linspace(time.t0,time.t0+time.simDt * time.stepNum,time.stepNum+1);
    
    %% parameter related to error
    % 初期時計誤差
    error.initialClock = 10   * randn; %初期時計誤差(秒)
    error.randomClock        = 0.01;         %ランダム時計誤差
    % 初期探査機軌道誤差[km]. 
    error.scPos0 = randn(3,1) *  1000;
    % 適当に0.1km/s程度の誤差とする
    error.scVel0 = randn(3,1) *  0.1;
    % ダイナミクスの不確定性の標準偏差(探査機)
    error.dynamics = 1e-10;
    
    %% ground station
    gs.lat  = 36.1325063*cspice_rpd();
    gs.lon =138.3607113*cspice_rpd();
    gs.alt  = 1.456;
    %% laser form gs to sc
    % 探索範囲(rad)
    gs.searchArea     = 5e-6; %±5μrad(だいぶ余裕を持たせている)
    gs.searchStep     = 1e-7; %探索時の1stepあたりの間隔(rad)
    gs.searchTimeStep = 2e-2; %適当．SOTAの資料にこれに相当するかは分からないが20msの記述あり
    % 探索1回にかかる時間
    time.obs = (2 * gs.searchArea/gs.searchStep)^2 * gs.searchTimeStep;
    % 探索1回にかかる時間がシミュレーションの何stepに相当するか
    time.obsStep = round(time.obs/time.simDt);
    % 地上局のレーザーに関するパラメータ
    gs.peakPower     = 370*10^3;
    gs.meanPower     = 2.5*10^3;
    gs.tEff  = 0.7;
    gs.atmosphereEff = 0.73;
    gs.aperture      = 1;
    gs.gamma         = 0.2; % telescope obscuration ratio
    gs.wavelength    = 1064 * 10^-9;
    gs.alpha         = 1.12 - 1.30 *  gs.gamma^2 + 2.12 * gs.gamma^4;
    gs.tAntGain      = (pi * gs.aperture/ gs.wavelength)^2 * 2/gs.alpha * (exp(-gs.alpha^2) - exp(-gs.alpha^2*gs.gamma^2))^2;
       
    %% 宇宙機の初期推定位置(とりあえず土星にしている)
    sc.state0 = cspice_spkezr('699', time.t0,'ECLIPJ2000', 'NONE', '10');
    % 光通信機
    sc.aperture = 0.22;
    sc.gamma    = 0.2;
    sc.rAntGain = (pi * gs.aperture/ gs.wavelength)^2 * 2/gs.alpha * (1- sc.gamma^2)^2;
    sc.rEff = 0.7;
    % QD センサ
    k = 1.381*1e-23;     %ボルツマン定数
    T = 298;             % 絶対温度
    Rsh = 50 * 10^6;     % 並列抵抗
    sc.qdBw = 2 * 1e7;   % 帯域幅
    sc.qdFov = 100*1e-6; % 視野角 
    sc.qdIj  = (4 * k * T * sc.qdBw/Rsh)^0.5; %熱雑音電流
    sc.qdS   = 0.68;     % 受光感度[A/W]
    sc.qdId  = 5*1e-10;   % 暗電流
    
    %% EKFの初期値, 観測誤差分散
    X_hat = [0;sc.state0];
    P = zeros(7);
    P(1,1) = 10^2;
    P(2,2) = 1000^2;
    P(3,3) = 1000^2;
    P(4,4) = 1000^2;
    P(5,5) = 1e-1^2;
    P(6,6) = 1e-1^2;
    P(7,7) = 1e-1^2;
    P_list = zeros(7,7,length(time.list));
    P_list(:,:,1) = P;
    
    R = [1e-6,0,0;      %観測誤差分散
         0,1e-12,0;
         0,0,1e-12];
    
    
    
    
   

end