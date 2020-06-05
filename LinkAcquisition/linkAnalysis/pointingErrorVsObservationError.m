% 地上局の指向誤差と探査機の観測誤差の関係を求める
%% constant value
%     constant.sunMu         = SSD.GM(10);      % km^3/s^2 gravity constant of the sun
%     constant.earthMu       = SSD.GM(399);     % km^3/s^2 gravity constant of the 
%     constant.saturnMu      = SSD.GM(799);     %[km/s^2]
%     constant.lightSpeed    = 299792.458;      % [km/s] light speed 
%     constant.earthRadius   = SSD.r(399);      % [km] earth radius
%     constant.earthRotation = 2*pi/24/60/60;   % [rad/s] earth rotation speed
%     constant.eathAxis      = 23.4/180*pi;     % inclination of the earth axis
%     constant.earthF = 1/298.257223560; % https://topex.ucsd.edu/geodynamics/14gravity1_2.pdf
    constant.elementaryCharge = 1.602 * 1e-19;



tPointingError = 5e-7:5e-8:2e-6;
rPointingError = zeros(1,length(tPointingError));

gs.peakPower     = 370*10^3;
gs.meanPower     = 2.5*10^3;
gs.tEff  = 0.7;
gs.atmosphereEff = 0.73;
gs.aperture      = 1;
gs.gamma         = 0.2; % telescope obscuration ratio
gs.wavelength    = 1064 * 10^-9;
gs.alpha         = 1.12 - 1.30 *  gs.gamma^2 + 2.12 * gs.gamma^4;
gs.tAntGain      = (pi * gs.aperture/ gs.wavelength)^2 * 2/gs.alpha * (exp(-gs.alpha^2) - exp(-gs.alpha^2*gs.gamma^2))^2;
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

Ls = (gs.wavelength/(4 * pi * (1.496e+8 * 1e3)))^2;

for i= 1:length(tPointingError)
    pError = tPointingError(i);
    Lp = calcPointingLoss(pError,gs);
    recPower = gs.peakPower * gs.tAntGain * gs.tEff * Lp * Ls * gs.atmosphereEff *  sc.rAntGain * sc.rEff; 
    rPointingError(i) = (sc.qdIj^2 ...
      + 2 * constant.elementaryCharge * sc.qdId * sc.qdBw ...
      + 2 * constant.elementaryCharge * recPower * sc.qdBw )^0.5 ...
      / (recPower * sc.qdS) * sc.qdFov;
  
end

semilogy(tPointingError,rPointingError)
xlabel("ground station transmit pointing error[rad]")
ylabel("gspacecraft observed pointing error[rad]")