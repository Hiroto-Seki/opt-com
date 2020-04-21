classdef groundState < handle
    properties
          dataType                %��ԗʂ̎��(����l, �^�l)
          t                       %����
          pos                     %position realtive to the central body [km]
          vel                     %velocity realtive to the central body [km]
          pos0
          vel0
    end
    methods 
        function obj = groundState(dataType) %�R���X�g���N�^ 
            obj.dataType = dataType;               
        end              
        getTrajectoryEarthRotation(obj,time,constant)
    end
    methods(Static)
        [gsPos,gsVel] = earthRotation(pos0, t, constant)
        [xv] = calcEarthRotation(xv, constant)
    end
end