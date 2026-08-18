function Gtx_dB = beamPattern(sat_lat, sat_lon, sat_alt,user_lat,user_lon,user_alt,cell_lat,cell_lon,cell_alt)

% Construct the circular-aperture antenna gain pattern.
c = 3e8;
f = 2e9;
lamda = c/f;

Nant=32;
Gt0=Nant^2;

k=2*pi*f/c;
a=1;

theta = -90:0.15:90;
J1 = besselj(1,k*a*sind(theta));

G = zeros(size(theta));

G(theta == 0) = 1;
G(theta ~= 0) = 4 * (abs(J1(theta ~= 0) ./ (k * a * sind(theta(theta ~= 0)))).^2);
G_lin=Gt0.*G;

% Determine the off-axis angle between the user and serving-cell directions.
wgs84 = wgs84Ellipsoid('kilometer');

[satX, satY, satZ] = geodetic2ecef(wgs84, sat_lat, sat_lon, sat_alt,"degrees");
[ueX, ueY, ueZ] = geodetic2ecef(wgs84, user_lat, user_lon, user_alt,"degrees");
[cellX, cellY, cellZ] = geodetic2ecef(wgs84, cell_lat, cell_lon, cell_alt,"degrees");

vectorSatToUE = [ueX - satX, ueY - satY, ueZ - satZ];
vectorSatToCell=[cellX-satX, cellY-satY, cellZ-satZ];

den = norm(vectorSatToUE) * norm(vectorSatToCell);
if den < 1e-12
    angle_AOB = 0;
else
    cos_AOB = dot(vectorSatToUE, vectorSatToCell) / den;
    cos_AOB = min(1, max(-1, cos_AOB));
    angle_AOB = acosd(cos_AOB);
end

% Interpolate the linear gain and return it in decibels.
Gtx = interp1(theta, G_lin, angle_AOB, 'makima', 'extrap');
Gtx = max(Gtx, realmin);

Gtx_dB = 10 * log10(Gtx);

end
