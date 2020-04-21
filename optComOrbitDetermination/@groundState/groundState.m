classdef groundState < handle
    properties
          dataType                %状態量の種類(推定値, 真値)
          t                       %時刻
          pos                     %position realtive to the central body [km]
          vel                     %velocity realtive to the central body [km]
          pos0
          vel0
    end
    methods 
        function obj = groundState(dataType) %コンストラクタ 
            obj.dataType = dataType;               
        end              
        getTrajectoryEarthRotation(obj,time,constant)
    end
    methods(Static)
        [gsPos,gsVel] = earthRotation(pos0, t, constant)
        [xv] = calcEarthRotation(xv, constant)
    end
end