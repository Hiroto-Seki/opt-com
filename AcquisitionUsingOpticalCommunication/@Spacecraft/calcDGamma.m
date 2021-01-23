function dGamma = calcDGamma(tau,delFdelX,dt)
    dGamma = expm(delFdelX * (dt-tau));
    dGamma = reshape(dGamma,[7^2,1]);
end