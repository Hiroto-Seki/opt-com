% get ephemeris of eath 
function getEphem(obj)
    obj.state = cspice_spkezr('399', obj.t,'ECLIPJ2000', 'NONE', '10');
end