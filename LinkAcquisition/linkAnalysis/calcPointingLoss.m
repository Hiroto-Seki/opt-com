% input 
% theta  : pointing error (rad)
% lambda : wavelength     ( m )
% D      :  diameter of transmitter
% gamma  : ratio of obscuration diameter to the primary aperture diameter
% (typical value = 0.2)

% output: Lp = pointing efficiency


function Lp = calcPointingLoss(theta,gs)  
    % 分子の積分
    intN = 0;
    % 分母の積分
    intD = 0;
    for u = gs.gamma^2:(1-gs.gamma^2)/20:1
        intN = intN + exp(-gs.alpha^2 * u) * besselj(0,pi*gs.aperture/gs.wavelength * theta * u^0.5);
        intD = intD + exp(-gs.alpha^2 * u);
    end
    Lp = (intN/intD)^2;
end