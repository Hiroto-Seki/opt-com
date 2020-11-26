% calculate laser transmit direction using gs estimates value
% 入力:
% i                : 地上局の時刻gsTrue.t(i)に対応する値を計算する
% gsTrue
% earth            : 地球
% time
% constant
% 出力
%  gsTrue.t_ut                %真値に一番近い方向を照射する時の時刻
%  gsTrue.direction_ut        %光照射方向
%  

function [gsTrue,earth] = search(i,gsTrue,earth,gs,time,constant) 
 % 探索方向は，慣性空間での角度．自身の速度は補正していない値．(実際に観測される値ではないが，地上局ならこれくらいできそう)
 % 推定値の方向の初期値(中心となる)
 relXvEst = gsTrue.opnEstTempState_ut - earth.state(:,i) - gsTrue.state(:,i); 
 relXvTrue = gsTrue.opnTrueTempState_ut - earth.state(:,i) - gsTrue.state(:,i);
 estAzm0 = atan2(relXvEst(2),relXvEst(1) );
 estElv0 = atan(relXvEst(3) /(relXvEst(1)^2 + relXvEst(2)^2)^0.5 );
 % 真値の方向の初期値
 trueAzm0  = atan2(relXvTrue(2),relXvTrue(1) );
 trueElv0 = atan(relXvTrue(3) /(relXvTrue(1)^2 + relXvTrue(2)^2)^0.5 );
 % 推定値の速度(一定とする)
 estAzmVel  = (relXvEst(5)*relXvEst(1) - relXvEst(2)*relXvEst(4))/ (relXvEst(1)^2 + relXvEst(2)^2);
 estElvVel = ( relXvEst(6)*(relXvEst(1)^2 + relXvEst(2)^2) - relXvEst(3)*(relXvEst(1)*relXvEst(4) + relXvEst(2)+ relXvEst(5)) )...
     /(relXvEst(1)^2 + relXvEst(2)^2 + relXvEst(3)^2) / (relXvEst(1)^2 + relXvEst(2)^2)^0.5;
 % 真値の速度(一定とする)
 trueAzmVel = (relXvTrue(5)*relXvTrue(1) - relXvTrue(2)*relXvTrue(4))/ (relXvTrue(1)^2 + relXvTrue(2)^2);
 trueElvVel = ( relXvTrue(6)*(relXvTrue(1)^2 + relXvTrue(2)^2) - relXvTrue(3)*(relXvTrue(1)*relXvTrue(4) + relXvTrue(2)+ relXvTrue(5)) )...
     /(relXvTrue(1)^2 + relXvTrue(2)^2 + relXvTrue(3)^2) / (relXvTrue(1)^2 + relXvTrue(2)^2)^0.5;
 
 % plot用
 estAzmList = zeros(1, round((2 * gs.searchArea/gs.searchStep)^2));
 estElvList = zeros(1, round((2 * gs.searchArea/gs.searchStep)^2));
 trueAzmList = zeros(1, round((2 * gs.searchArea/gs.searchStep)^2));
 trueElvList = zeros(1, round((2 * gs.searchArea/gs.searchStep)^2));
 pointAzmList = zeros(1, round((2 * gs.searchArea/gs.searchStep)^2));
 pointElvList = zeros(1, round((2 * gs.searchArea/gs.searchStep)^2));
 
 
 % 探索に関するパラメータ
 i_azm = 0;
 i_elv = 0;
 i_dir = 1; %次の探索方向　1=+azm, 2=+elv, 3=-azm, 4=-elv
 dt = 0;
 tempPointingError = 1;
 for j = 1:(2 * gs.searchArea/gs.searchStep)^2
     % 各時刻のAzmとElvを計算．(推定値真値ともに)
     estAzm = estAzm0 + dt * estAzmVel;
     estElv = estElv0 + dt * estElvVel;
     pointAzm = estAzm + i_azm * gs.searchStep;
     pointElv = estElv + i_elv * gs.searchStep;
     trueAzm = trueAzm0 + dt * trueAzmVel;
     trueElv = trueElv0 + dt * trueElvVel;
     % plot用に格納
     estAzmList(j) = estAzm;
     estElvList(j) = estElv;
     pointAzmList(j) = pointAzm;
     pointElvList(j) = pointElv;
     trueAzmList(j) = trueAzm;
     trueElvList(j) = trueElv;
     
     % 各時刻の角度誤差を計算する
     azmError = min(mod(pointAzm -trueAzm, 2*pi ),mod( - pointAzm + trueAzm, 2*pi));
     elvError = min(mod(pointElv -trueElv, 2*pi ),mod( - pointElv + trueElv, 2*pi));
     pointingError = (azmError^2 + elvError^2)^0.5;
     % 一番pointingErrorが小さくなったら更新する．
     if pointingError < tempPointingError
         tempPointingError = pointingError;
         gsTrue_t_ut   = gsTrue.t(i) + dt;
         gsTrue_direction_ut = [pointAzm;pointElv];
     end
     % 次時刻ごとに推定値から何step離れた点に照射するか計算
     if     mod(i_dir,4) == 1
         i_azm = i_azm + 1;
         if round(j - ((i_dir+1)/2)^2) ==0
             i_dir = i_dir +1;
         end
     elseif mod(i_dir,4) ==2
         i_elv = i_elv + 1;
          if round(j - (i_dir/2)^2 - i_dir/2) ==0
             i_dir = i_dir +1;
         end
     elseif mod(i_dir,4) ==3
         i_azm = i_azm -1;
         if round(j - ((i_dir+1)/2)^2) ==0
             i_dir = i_dir +1;
         end
     else 
         i_elv = i_elv -1;
         if round(j - (i_dir/2)^2 - i_dir/2) ==0
             i_dir = i_dir +1;
         end
     end
     dt = dt + gs.searchTimeStep;
 end
 


 %% 何番目の送信に対応した状態量なのかを記録する
 gsTrue.ut_counter = gsTrue.ut_counter + 1;
 gsTrue.t_ut(gsTrue.ut_counter) = gsTrue_t_ut;
 gsTrue.direction_ut(:,gsTrue.ut_counter) = gsTrue_direction_ut;
  % 送信時刻の地上局の位置速度，地球の位置,速度を求める.
 earth.state_ut(:,gsTrue.ut_counter)  = earth.calcStateAtT_cb(gsTrue.t_ut(gsTrue.ut_counter),time);
 gsTrue.state_ut(:,gsTrue.ut_counter) = gsTrue.calcStateAtT_gs(gsTrue.t_ut(gsTrue.ut_counter),time,constant);
 

 
end