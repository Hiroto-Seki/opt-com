% 宇宙機がuplinkから受信する内容を求める
% obj = scTrue
function obj = receiveUplink(obj,gsTrue,earth,constant,time) 
    %% gsTrue.ut_counter番目のuplinkを宇宙機が受信する時刻と状態量を求める
    [obj.t_ur(gsTrue.ut_counter),obj.state_ur(:,gsTrue.ut_counter)] ...
        = GroundStation.calcTarget(gsTrue.t_ut(gsTrue.ut_counter),gsTrue.state_ut(:,gsTrue.ut_counter),earth.state_ut(:,gsTrue.ut_counter),obj,time,constant);  
    
    %% 受信する内容
    obj.gsT_ur(gsTrue.ut_counter)               =  gsTrue.t_ut(gsTrue.ut_counter);              %uplinkに載っている，送信時刻
    obj.eState_ur(:,gsTrue.ut_counter)          =  earth.state_ut(:,gsTrue.ut_counter);            %uplinkに載っている，送信時刻の地球の状態量   
    obj.gsState_ur(:,gsTrue.ut_counter)         =  gsTrue.state_ut(:,gsTrue.ut_counter);           %uplinkに載っている, 送信時刻の地上局の状態量
    obj.transDirection_ur(:,gsTrue.ut_counter)  =  gsTrue.direction_ut(:,gsTrue.ut_counter);       %uplinkに載っている，送信方向(誤差を持つ)
end