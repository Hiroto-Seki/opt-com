% get ephemeris of celestial body 
function getEphem(obj)
    if strcmp(obj.name,"Earth")
    obj.state = cspice_spkezr('399', obj.t,'ECLIPJ2000', 'NONE', '10'); %太陽を中心とした地球
    end
end