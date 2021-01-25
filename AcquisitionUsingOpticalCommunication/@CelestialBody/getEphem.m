% get ephemeris of celestial body 
function getEphem(obj,time)
    if strcmp(obj.name,"Earth")
    obj.state = cspice_spkezr('399', time.t0Ephemeris + obj.t,'ECLIPJ2000', 'NONE', '10'); %太陽を中心とした地球
    elseif strcmp(obj.name,"Saturn")
        obj.state = cspice_spkezr('699', time.t0Ephemeris + obj.t,'ECLIPJ2000', 'NONE', '10'); %太陽を中心とした地球
    end
end 