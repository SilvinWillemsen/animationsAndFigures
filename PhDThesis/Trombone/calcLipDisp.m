%% Obtain deltaP
a1 = 2 / k + omega0^2 * k + sig + g^2 * k / (2 * M);
a2 = Sr / M;
a3 = 2/k * 1/k * (y - yPrev) - omega0^2 * yPrev + g / M * psiPrev;
b1 = SHalf(1) * vNext(1) + h * SBar(1) / (rho * c^2 * k) * (Pm  - p(1));
b2 = h * SBar(1) / (rho * c^2 * k);
c1 = w * subplus(y + H0) * sqrt(2 / rho);
c2 = b2 + a2 * Sr / a1;
c3 = b1 - a3 * Sr / a1;

deltaP = sign(c3) * ((-c1 + sqrt(c1^2 + 4 * c2 * abs(c3)))/ (2 * c2))^2;

%% Update lip scheme
gammaR = g * k^2 / (2 * M);
alpha = 2 + omega0^2 * k^2 + sig * k + g * gammaR;
beta = sig * k - 2 - omega0^2 * k^2 + g * gammaR;
xi = 2 * Sr * k^2 / M;

yNext(n) = 4 / alpha * y + beta / alpha * yPrev + xi / alpha * deltaP + 4 * gammaR * psiPrev / alpha;
