function ret = scalingFactor(a, q)
%Return a row vector representing phi(q), the scaling factor
%phi(q) = F(q)^2 / q * absorption
%D is the D-spacing at which F(q) was obtained.
%a is a vector for the coefficients of each cosine term
%q is a vector whose elements are the values of q for each phi(q)
global largeNum;
if iscolumn(a) == 1
    a = transpose(a);
end
if isrow(q) == 1
    q = transpose(q);
end
D = 62.9;
constraint = 0;
%cosineMatrix: each column corresponds to a term in the sum in F(q)
cosineMatrix = zeros(length(q),length(a));
for k = 1 : length(a)
    cosineMatrix(:,k) = a(k) * (-1)^k * 2 * sin(q*D/2) .* (-1 ./ q + D^2 * q ./ (q.^2 * D^2 - 4 * pi^2 * k^2));
end
%Sum all the cosine terms to get F(q)
ret = sum(cosineMatrix,2);
%Square F(q)
ret = ret .* ret;
%Apply Lorentz polarization factor appropriate for oriented sample
ret = ret ./ q;
%Apply absorption factor

%Return scaling factor, phi(q), as a row vector
ret = transpose(ret);
%Check whether the derived electron density profile is reasonable. If not, the last element in ret
%will be changed to a large value. This element gets compared to a number, zero, and the difference
%contributes to the chisquare value computed by lsqcurvefit.
zstep = 0.5;
z = [-D/2:zstep:D/2];
rho = edp(a, z);
temp = isreasonable(z, rho);
if temp == false
    constraint = largeNum;
end
%ret(length(ret)) = constraint;

end

