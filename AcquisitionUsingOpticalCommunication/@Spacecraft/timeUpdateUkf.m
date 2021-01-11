% UKFでの，シグマ点列，共分散行列を時間更新する
% 入力 
% x_sp: 現時点でのシグマ点列
% ukf: ukfに関するパラメータ
% Dt: 更新する時間
% dt: タイムステップ(Dt>dt)の時，この間隔で伝搬する
function [X_new, P_new, x_sp_new] = timeUpdateUkf(x_sp,constant, ukf,Dt,dt, error)
% 初期化
X_new = zeros(size(x_sp,1),1);
P_new = zeros(size(x_sp,1));
x_sp_new = zeros(size(x_sp,1),size(x_sp,2));

% Qの計算
Q = [zeros(4,7);zeros(3,4),eye(3)*(error.dynamics * Dt)^2];

% シグマ点列の時間更新
for i = 1:size(x_sp,2)
   % 初期化
   xvsc = x_sp(2:7,i);
   xvsc = Spacecraft.timeUpdate_sc(xvsc,constant.sunMu, Dt, dt);
   % シグマ点列の代入
   x_sp_new(1,i)   = x_sp(1,i);
   x_sp_new(2:7,i) = xvsc;
end
% 更新したシグマ点列から状態ベクトルの平均を求める
for j = 1:size(x_sp,2)
    if j==1
        X_new = X_new + ukf.w0_m * x_sp_new(:,j);
    else
        X_new = X_new + ukf.wi_m * x_sp_new(:,j);
    end  
end
% 誤差共分散行列を求める
for k = 1:size(x_sp,2)
    if k==1
        P_new = P_new + ukf.w0_c * (x_sp_new(:,k) - X_new) * (x_sp_new(:,k) - X_new).';
    else
        P_new = P_new + ukf.wi_c * (x_sp_new(:,k) - X_new) * (x_sp_new(:,k) - X_new).';
    end  
end


P_new = P_new + Q;

end