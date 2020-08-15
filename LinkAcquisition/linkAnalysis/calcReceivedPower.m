% 起動決定精度[m]と受信信号レベル[W]の比較
% uplinkを想定

clear; 

% 定数
P_l_mean   = 2500;           % 送信電力[W] 
P_l_peak   = 370*10^3;       % ピーク送信電力[W]
lambda     = 1064 * 10^(-9); % 波長[m]
D_t        = 1;              % 送信口径[m] 
D_r        = 0.22;           % 受信口径[m] 
gamma      = 0.2;            % obscuration ratio[]
eta_t      = 0.7;            % 送信効率 []
eta_r      = 0.7;            % 受信効率 []
l_a        = 0.73;           % 大気損失 []
eta_lambda = 1;              % バンドフィルター損失 []
eta_D      = 1;              % 受信指向誤差損失 []
L          = 1.5*10^12;      % 伝搬距離[m] 

reqPower   = 10^(-13);       % 受信信号レベルの下限値

alpha_t    = 1.12 -1.30*gamma^2 + 2.12*gamma^4;
% 送信アンテナ利得
g_t        = (pi*D_t/lambda)^2 *2/alpha_t * (exp(-alpha_t^2) - exp(-alpha_t^2*gamma^2))^2;
% 受信アンテナ利得
g_r        = (pi*D_r/lambda)^2 *2/alpha_t * (1-gamma^2);
% 自由空間損失
l_s = (lambda/(4*pi*L))^2;

% 探査機の位置精度(m)
posErrorList = 0:10^4:5*10^6;
P_r_meanList = zeros(1, length(posErrorList));
P_r_peakList = zeros(1, length(posErrorList));

for i = 1:length(posErrorList)
    theta = posErrorList(i)/L;
    eta_tp = calcPointingLoss(theta,gamma,alpha_t,D_t,lambda);
    P_r_mean = P_l_mean * g_t * eta_t * eta_tp* l_s * l_a * g_r * eta_r * eta_lambda * eta_D;
    P_r_peak = P_l_peak * g_t * eta_t * eta_tp* l_s * l_a * g_r * eta_r * eta_lambda * eta_D;
    P_r_meanList(1,i) = P_r_mean;
    P_r_peakList(1,i) = P_r_peak;
end

hold on
plot(posErrorList/10^3, 10*log10(P_r_peakList*10^3))
plot(posErrorList/10^3, 10*log10(P_r_meanList*10^3))
yline(10*log10(reqPower*10^3*10^0.3))
xlabel('position error[km]')
ylabel('received power[dBm]')
legend('received peak power', 'received mean power', 'required power+3dB')
hold off






