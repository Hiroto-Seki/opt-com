% STTによる観測と，QDセンサーによる観測から，その時刻の姿勢を推定する．(最小分散推定)
% 宇宙機の軌道とカップリングがあるが，ここでは，軌道の推定値の更新は行わない
% ---- 入力 ----
% stt: sttによる観測値(roll,pitch,yaw)→これをノミナルの状態量に使う
% sttError: 観測誤差分散
% qd : qdによる観測値(qd_x, qd_y)
% directionEstI :
% 宇宙機の軌道の推定値を用いて計算した，慣性空間上での見かけのuplink受信方向．qdセンサーの出力を慣性空間上の値に変換するのに使う
% scfl: 光学系の焦点距離. 単位がm
%---- 出力 ----
% roll, pitch,yaw: 推定した宇宙機の姿勢

%---- 中間パラメーター -----
% X      : 最終的な推定値
% X_star : アプリオリな推定値
% x_hat  : X - X_star
% Y      : 観測値
% Y_star : アプリオリな推定値に対応する観測値
% y      : Y - Y_star
% H      : 観測方程式を状態量(=アプリオリな推定値)で微分したもの
% R      : 観測誤差共分散行列
% P      : 誤差共分散行列

function [roll, pitch, yaw, P] = attDetermination(stt,sttError,qd,qdError,directionEstI,scfl)
  % 観測値をまとめる
  Y = [stt;qd];
  % アプリオリな推定値を与える
  X_star = stt;
  %  観測誤差共分散行列を与える
  R = [sttError^2 * eye(3),         zeros(3,2);
                zeros(2,3), qdError^2 * eye(2)];
            
  % アプリオリな状態量での回転行列を求める．
  rotation_star = Spacecraft.rotation(-X_star(1),-X_star(2),-X_star(3),2);
  % アプリオリな状態量での光学系座標系での見かけのuplink方向を求める
  directionFlt_star = rotation_star * directionEstI;
  % アプリオリな状態量でのqdセンサーの観測値を求める
  qd_star = - scfl/(directionFlt_star(3)*1e3) * [1, 0, 0; 0, 1, 0] * (directionFlt_star * 1e3);  % 単位がmであることに注意
  % Y_starを求める
  Y_star = [stt;qd_star];
  % yを求める
  y = Y - Y_star;
  
  %% Hを求める
  % 回転行列(rotation_star)の各要素
  R_x = [               1,               0,               0;
                        0,  cos(X_star(1)),  sin(X_star(1));
                        0, -sin(X_star(1)),  cos(X_star(1))];
  R_y = [  cos(X_star(2)),               0, -sin(X_star(2));
                        0,               1,               0;
           sin(X_star(2)),               0,  cos(X_star(2))];
  R_z = [  cos(X_star(3)),  sin(X_star(3)),               0;
          -sin(X_star(3)),  cos(X_star(3)),               0;
                        0,               0,               1];
  %回転行列の各要素の微分
  difR_x = [               0,               0,               0;
                           0, -sin(X_star(1)),  cos(X_star(1));
                           0, -cos(X_star(1)), -sin(X_star(1))];
  difR_y = [ -sin(X_star(2)),               0, -cos(X_star(2));
                           0,               1,               0;
              cos(X_star(2)),               0, -sin(X_star(2))];
  difR_z = [ -sin(X_star(3)),  cos(X_star(3)),               0;
             -cos(X_star(3)), -sin(X_star(3)),               0;
                           0,               0,               1];
  % 回転行列の微分
  delRdelPhi   = difR_x * R_y * R_z;
  delRdelTheta = R_x * difR_y * R_z;
  delRdelPsi   = R_x * R_y * difR_z;
  % qdセンサーの観測方程式をphi(x軸周りの回転)で微分する %単位がmであることに注意
  delQdDelPhi    = -1/(directionFlt_star(3)*1e3)^2 *  [0,0,1] * delRdelPhi * (directionEstI * 1e3) *...
      (-scfl) * [1, 0, 0; 0, 1, 0] * (directionFlt_star*1e3)...
      - scfl/(directionFlt_star(3)*1e3) * [1, 0, 0; 0, 1, 0] * delRdelPhi * (directionEstI * 1e3);
  delQdDelTheta  = -1/(directionFlt_star(3)*1e3)^2 *  [0,0,1] * delRdelTheta * (directionEstI * 1e3) *...
      (-scfl) * [1, 0, 0; 0, 1, 0] * (directionFlt_star * 1e3)...
      - scfl/(directionFlt_star(3)*1e3) * [1, 0, 0; 0, 1, 0] * delRdelTheta * (directionEstI*1e3);
  delQdDelPsi  = -1/(directionFlt_star(3)*1e3)^2 *  [0,0,1] * delRdelPsi * (directionEstI*1e3) *...
      (-scfl) * [1, 0, 0; 0, 1, 0] * (directionFlt_star*1e3)...
      - scfl/(directionFlt_star(3)*1e3) * [1, 0, 0; 0, 1, 0] * delRdelPsi * (directionEstI*1e3);
  % Hを求める
  H = [eye(3);delQdDelPhi,delQdDelTheta,delQdDelPsi];
  
  % Pを求める
  P = inv(H.' * inv(R) * H);
  
  %% x_hatを求める
  x_hat = P * H.' * inv(R) * y;
  
  X = X_star + x_hat;
  
  roll  = X(1);
  pitch = X(2);
  yaw   = X(3);
  

end