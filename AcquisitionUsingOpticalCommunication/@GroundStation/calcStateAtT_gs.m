% ある時刻tの地上局の状態量を計算する
function  xvAtT = calcStateAtT_gs(obj,t,time,constant)
% 最終時刻までのリスト内に収まっている時
    if t < obj.t(length(obj.t))
        closeTimeIndex = floor( (t - obj.t(1) )/ time.simDt ) + 1;
        closeTimeOffset = t - obj.t(closeTimeIndex);      % 一番近い時刻から伝搬しなければいけない時間
        xvg = obj.state(:,closeTimeIndex);
        xvAtT = GroundStation.earthRotation(xvg(1:3), closeTimeOffset, constant);
     % 最終時刻までのリストに収まっていない時
    else
        restTime = t -  obj.t(length(obj.t)); 
        xvg =  obj.state(:,length(obj.t));
        xvAtT = GroundStation.earthRotation(xvg(1:3), restTime, constant);
    end
end