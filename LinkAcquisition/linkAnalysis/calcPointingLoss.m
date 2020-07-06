% input 
% theta  : pointing error (rad)
% lambda : wavelength     ( m )
% D      :  diameter of transmitter
% gamma  : ratio of obscuration diameter to the primary aperture diameter
% (typical value = 0.2)

% output: Lp = pointing efficiency


function Lp = calcPointingLoss(theta,gamma,alpha,aperture,wavelength)  
    % 分子の積分
    intN = 0;
    % 分母の積分
    intD = 0;
    for u = gamma^2:(1-gamma^2)/20:1
        intN = intN + exp(-alpha^2 * u) * besselj(0,pi*aperture/wavelength * theta * u^0.5);
        intD = intD + exp(-alpha^2 * u);
    end
    Lp = (intN/intD)^2;
end