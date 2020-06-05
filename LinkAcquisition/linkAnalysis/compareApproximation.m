% pointing loss のグラフとその近似式のグラフを作る
% 近似式はgamma=0.2のもの
% うまく近似できる範囲を確認するためのもの

gamma = 0.2;
D = 1;
lambda = 10^-6;
thetaList = 0:0.05*lambda/D:lambda/D*3;
LpList = zeros(1,length(thetaList));

% 近似値の計算に用いる
f0 = 0.555645;
f2 = -0.120457;
f4 = 0.0542465;
f6 = -0.00317773;

for i = 1:length(thetaList)
    theta = thetaList(i);
    Lp = calcPointingLoss(theta, lambda, D, gamma);
    LpList(1,i) = Lp;
    LpAprx = 1/f0^2*( f0 + f2/2 * (pi*D/lambda * theta)^2 ...
        + f4/(2*3*4) * (pi*D/lambda * theta)^4 + f6/(2*3*4*5*6) * (pi*D/lambda * theta)^6  )^2;
    LpAprxList(1,i) = LpAprx;
end

hold on
plot(thetaList, 10*log10(LpList))
plot(thetaList, 10*log10(LpAprxList))
legend('exact','approximation')
xlabel('pointing error[rad]')
ylabel('pointing loss[dB]')
hold off
