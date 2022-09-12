
clear all

z = 3 * rand(1) + 3i * randn(1)

Gamma = (z-1) / (z+1)
a = abs(Gamma)
th = angle(Gamma)

zCircleCenter = (1+a^2) / (1-a^2)
zCircleRadius = 2*a / (1-a^2)
zCircleAngle  = 2 * atan(sin(th) ./ (cos(th) - a)) - th

zSynth = zCircleCenter + zCircleRadius * exp(1i*zCircleAngle)

error = zSynth - z
