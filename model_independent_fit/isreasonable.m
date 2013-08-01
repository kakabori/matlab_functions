function ret = isreasonable(z, rho)
%Return ture if the electron density profile is reasonable. Return false otherwise.
%   Detailed explanation goes here

ret = true;
% Look for the first peak position from the left (negative z side)
[value index] = max(rho);
% If the maximum has a negative value, unreasonable electron density profile is obtained. For
% loop continues on to the next iteration and the current fit result does not get saved.
if value <= 0
    ret = false;
    return;
end
% The peak position should correspond to a phosphate group, which had better be between -25 and
% -14 Angstroms. These values should be optimized for a given data set.
if z(index) < -25 && z(index) > - 14
    ret = false;
    return;  
end
% Look for the terminal trough.
[value index] = min(rho);
% If the minimum is not negative, the fit is unreasonable.
if value >= 0
    ret = false;
    return;
end
% If the trough is not at the bilayer center, the fit is unreasonable.
if z(index) ~= 0
    ret = false;
    return
end

end

