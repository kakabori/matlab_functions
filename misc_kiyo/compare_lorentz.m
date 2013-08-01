function [ms1 ms2 ms3] = compareLorentz(gamma)

x = (-50*gamma):(gamma/10):(50*gamma);

figure;
hold on;
plot(x, yLorentz(gamma, x, 1*gamma), 'b');
plot(x, gamma2Lorentz(gamma, x, 1*gamma), 'r');
plot(x, yLorentz(gamma, x, 5*gamma), 'b');
plot(x, gamma2Lorentz(gamma, x, 5*gamma), 'r');
plot(x, yLorentz(gamma, x, 10*gamma), 'b');
plot(x, gamma2Lorentz(gamma, x, 10*gamma), 'r');

temp = (yLorentz(gamma, x, 1*gamma) - gamma2Lorentz(gamma, x, 1*gamma)) .^ 2;
ms1 = sum(temp);
temp = (yLorentz(gamma, x, 5*gamma) - gamma2Lorentz(gamma, x, 5*gamma)) .^ 2;
ms2 = sum(temp);
temp = (yLorentz(gamma, x, 10*gamma) - gamma2Lorentz(gamma, x, 10*gamma)) .^ 2;
ms3 = sum(temp);


end

function ret = yLorentz(gamma, x, y)
ret = (gamma / pi) ./ (x.^2 + y^2 + gamma^2);
end