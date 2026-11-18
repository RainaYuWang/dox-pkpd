function [CF, CB] = fig3_pk(t, dose, wash, kBF)
% Exact solution for the fixed linear PK system and a rectangular CE input.
kFE = 1.25; kFB = 0.21;
input = (3.5e-6 * 250 / 0.004629) * dose;
total = kFE + kFB + kBF;
gap = sqrt((kFE+kFB-kBF)^2 + 4*kFB*kBF);
fast = (total + gap) / 2;
slow = kFE*kBF / fast;
u = min(t, wash);
s = max(0, t-wash);
bound = input*kFB/gap .* (integralExp(slow,u) - integralExp(fast,u));
slope = input*kFB/gap .* exp(-slow*u) .* (-expm1(-gap*u));
% Propagate the washout state. This includes both return flux and rebinding.
factor = exp(-slow*s);
mix = -expm1(-gap*s) / gap;
CB = factor .* (bound .* (1 + slow*mix) + slope.*mix);
dCB = factor .* (slope .* exp(-gap*s) - (slow*(fast*bound+slope)).*mix);
CF = (dCB + kBF*CB) / kFB;
end

function v = integralExp(rate,t)
if rate == 0
    v = t;
else
    v = -expm1(-rate*t) / rate;
end
end
