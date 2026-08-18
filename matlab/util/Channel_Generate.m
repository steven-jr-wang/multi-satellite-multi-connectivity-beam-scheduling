function h = Channel_Generate(sat_lat, sat_lon, sat_alt,user_lat,user_lon,user_alt,cell_lat,cell_lon,cell_alt)

% Compute link distance and elevation angle in ECEF coordinates.
wgs84 = wgs84Ellipsoid('kilometer');

[satX, satY, satZ] = geodetic2ecef(wgs84, sat_lat, sat_lon, sat_alt,"degrees");
[ueX, ueY, ueZ] = geodetic2ecef(wgs84, user_lat, user_lon, user_alt,"degrees");

vectorSatToUE = [ueX - satX, ueY - satY, ueZ - satZ];
distanceSatToUE = norm(vectorSatToUE);

distanceSatToCenter = norm([satX, satY, satZ]);

distanceUEToCenter = norm([ueX, ueY, ueZ]);

input_elevation_angle=asind(-((distanceSatToUE^2 + distanceUEToCenter^2 - distanceSatToCenter^2) / (2 * distanceSatToUE * distanceUEToCenter)));

% Randomly select the suburban or rural channel scenario.
scenario_index = randi([1, 2]);

elevation_angles = [10, 20, 30, 40, 50, 60, 70, 80, 90];
suburban_rural_los = [78.2, 91.9, 91.9, 92.9, 93.5, 94.0, 94.9, 95.2, 99.8];

suburban_rural_prob = interp1(elevation_angles, suburban_rural_los, input_elevation_angle, 'spline');

los_rand=rand;
kai=(los_rand<=0.01*suburban_rural_prob);

% Generate free-space, shadow, and clutter losses.
f=2e9;
A_FS=32.45+20*log10((f/(1e6)))+20*log10(distanceSatToUE);

if kai==1

    LOS_SF = [1.79, 1.14, 1.14, 0.92, 1.42, 1.56, 0.85, 0.72, 0.72];

    A_SF = interp1(elevation_angles, LOS_SF, input_elevation_angle, 'spline');

    A_total=A_FS+A_SF;

else

    NLOS_SF = [8.93, 9.08, 8.78, 10.25, 10.56, 10.74, 10.17, 11.52, 11.52];
    S_band_CL = [19.52, 18.17, 18.42, 18.28, 18.63, 17.68, 16.50, 16.30, 16.30];
    A_SF = interp1(elevation_angles, NLOS_SF, input_elevation_angle, 'spline');
    A_CL = interp1(elevation_angles, S_band_CL, input_elevation_angle, 'spline');

    A_total=A_FS+A_SF+A_CL;
end

Gtx=beamPattern(sat_lat, sat_lon, sat_alt,user_lat,user_lon,user_alt,cell_lat,cell_lon,cell_alt);
Grx=0;

ksi=sqrt(10^((Gtx+Grx-A_total)/10));

% Select scenario- and LOS-dependent small-scale fading parameters.
if scenario_index==1
    if kai==1

        DS_mu = [-8.16, -8.56, -8.72, -8.71, -8.72, -8.66, -8.38, -8.34, -8.34];
        DS_sigma = [0.99, 0.96, 0.79, 0.81, 1.12, 1.23, 0.55, 0.63, 0.63];

        K_mu_dB_values = [11.40,19.45,20.80,21.20,21.60,19.75,12.00,12.85,12.85];
        K_sigma_dB_values=[6.26,10.32,16.34,15.63,14.22,14.19,5.70,9.91,9.91];

        r_tao_values=[2.20,3.36,3.50,2.81,2.39,2.73,2.07,2.04,2.04];

        if input_elevation_angle<=60
            nbrOfClusters=3;
        else
            nbrOfClusters=2;
        end
        nbrOfRayperCluster=20;
    else

        DS_mu = [-7.91, -8.39, -8.69, -8.59, -8.64, -8.74, -8.98, -9.28, -9.28];
        DS_sigma = [1.42, 1.46, 1.46, 1.21, 1.18, 1.13, 1.37, 1.50, 1.50];

        r_tao_values=[2.28, 2.33, 2.43, 2.26, 2.71, 2.10, 2.19, 2.06, 2.06];

        if input_elevation_angle<=50
            nbrOfClusters=4;
        else
            nbrOfClusters=3;
        end
        nbrOfRayperCluster=20;

    end
elseif scenario_index==2
    if kai==1

        DS_mu = [-9.55, -8.68, -8.46, -8.36, -8.29, -8.26, -8.22, -8.2, -8.19];
        DS_sigma = [0.66, 0.44, 0.28, 0.19, 0.14, 0.1, 0.1, 0.05, 0.06];

        K_mu_dB_values = [24.72, 12.31, 8.05, 6.21, 5.04, 4.42, 3.92, 3.65, 3.59];
        K_sigma_dB_values=[5.07, 5.75, 5.46, 5.23, 3.95, 3.75, 2.56, 1.77, 1.77];

        r_tao_values=[3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8];

        nbrOfClusters=2;
        nbrOfRayperCluster=20;
    else

        DS_mu = [-9.01, -8.37, -8.05, -7.92, -7.92, -7.96, -7.91, -7.79, -7.74];
        DS_sigma = [1.59, 0.95, 0.92, 0.92, 0.87, 0.87, 0.82, 0.86, 0.81];

        r_tao_values=[1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7];

        if input_elevation_angle<=20
            nbrOfClusters=3;
        else
            nbrOfClusters=2;
        end
        nbrOfRayperCluster=20;

    end

else 
    error('scenario_index must be 1 or 2.');
end

% Draw delay spread and Rician K-factor values.
DS_mu_interp = interp1(elevation_angles, DS_mu, input_elevation_angle, 'spline');
DS_sigma_interp = interp1(elevation_angles, DS_sigma, input_elevation_angle, 'spline');

log_DS = DS_mu_interp + DS_sigma_interp * randn();

DS = 10^log_DS;

if kai==1

    K_mu_dB_interp = interp1(elevation_angles, K_mu_dB_values, input_elevation_angle, 'spline');
    K_sigma_dB_interp = interp1(elevation_angles, K_sigma_dB_values, input_elevation_angle, 'spline');
    K_dB = K_mu_dB_interp + K_sigma_dB_interp * randn();

    K = 10^(K_dB / 10);

else

    K=0;

end

r_tao_interp = interp1(elevation_angles, r_tao_values, input_elevation_angle, 'spline');

% Generate cluster delays and normalized cluster powers.
tao_tmp=-r_tao_interp*repmat(DS,1,nbrOfClusters).*log(rand(1,nbrOfClusters));
tao_tmp=sort(tao_tmp-repmat(min(tao_tmp),1,nbrOfClusters));

zeta_values=3;
zeta=10^(zeta_values/10);
P_tmp = exp(-(r_tao_interp-1) * tao_tmp ./ (r_tao_interp * repmat(DS', 1, nbrOfClusters))) .* 10.^(-randn(1, nbrOfClusters) * (zeta^2)/10);

P = P_tmp ./ repmat(sum(P_tmp), 1, nbrOfClusters);

P = P ./ repmat((K + 1)', 1, nbrOfClusters);
P1 = [P,(K ./ (K + 1))'] ;

ind_abandon_NLOS = find(P1 ./ max(P1) <= 10^(-2.5));
P1(ind_abandon_NLOS) = 0;

P=P1(1:end-1);

ro=repmat(P,1,nbrOfRayperCluster);
P_ray=ro./nbrOfRayperCluster;

% Combine diffuse rays with the LOS component and large-scale gain.
realPart = sqrt(1/2) .* randn(1, nbrOfClusters*nbrOfRayperCluster);
imagPart = sqrt(1/2) .* randn(1, nbrOfClusters*nbrOfRayperCluster);
beta = realPart + 1j * imagPart;

small_scale_gain=sqrt(P_ray).*beta;

omega=abs(sqrt(K ./ (K + 1))+sum(small_scale_gain));

h=ksi*omega;

end
