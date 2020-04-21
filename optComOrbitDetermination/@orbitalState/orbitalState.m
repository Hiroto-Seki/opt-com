classdef orbitalState < handle
    properties
          name                    %name
          dataType                %状態量の種類(推定値, 真値)
          t                       %時刻
          centralBody             %ダイナミクスの中心天体
%           a                       %semi-major axis           (km)
%           e                       %Eccentricity Magnitude
%           i                       %inclination[rad]
%           AN                      %ascending node            (rad)
%           AP                      %argument of periapses     (rad)
%           f                       %true anomaly angle (rad)
          mu                      %Gravitational Constant of central body
          pos                     %position realtive to the central body [km]
          vel                     %velocity realtive to the central body [km]
          pos0
          vel0
          state0
          clockErrorCorrection
    end
    methods 
        function obj = orbitalState(name,dataType,centralBody) %コンストラクタ 
            obj.name = name;
            obj.dataType = dataType;
            obj.centralBody = centralBody;
                     
        end              
        elem2rv(obj);
        getOrbitTwoBody(obj,time,constant)
        calcOrbitTwoBody_pertubation(obj,time,constant, error)
%         getEphemData(obj,time,error)
        
    end
    methods(Static)
        xv = twobody(xv, mu,a_pertubation)
        A  = delFdelX(xv,mu)
        
    end
end