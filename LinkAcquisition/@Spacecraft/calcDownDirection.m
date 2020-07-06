% 概要: time.list分の探査機の位置を伝搬して求めていく
function calcDownDirection(obj,t,gs,e,scAtT,time,constant,j)
% 時刻tでの地上局位置，地球位置と，探査機の軌道から，探査機に光が届く時刻とその時の探査機の状態量を求める
   opnDown = Spacecraft.calcTarget(t,gs,e,scAtT,time,constant);
% 相対方向を求める．(相対速度込み)
  xve = opnDown.stateE;
  xvg = opnDown.stateGs;
  xvs = scAtT;
  ltd = opnDown.t - t;
  obj.azmDown(j)  =  atan2(xve(2) + xvg(2) - xvs(2) + xvs(5)*ltd...
                                , xve(1) + xvg(1) - xvs(1) + xvs(4)*ltd);
  obj.elvDown(j)  = atan( (xve(3) + xvg(3) - xvs(3) + xvs(6)*ltd)/...
                                   ((xve(2) + xvg(2) - xvs(2) + xvs(5) *ltd)^2 ...
                                  + (xve(1) + xvg(1) - xvs(1) + xvs(4) *ltd)^2)^0.5 );
                              
                              
  obj.tDown(j) = opnDown.t;
  obj.eDown(:,j) = xve;
  obj.gsDown(:,j) = xvg;

end
