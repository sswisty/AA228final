% Load the data

states = load("states.csv");
r_sat = states(:,1:3);
r_deb = states(:,7:9);

% Nominal Orbits
sat_nom = load("sat_nominal.csv");
deb_nom = load("deb_nominal.csv");

% Final orbits
sat_end = load("sat_end.csv");


figure; hold all
earthPlot;
plot3(sat_nom(:,1),sat_nom(:,2),sat_nom(:,3))
plot3(deb_nom(:,1),deb_nom(:,2),deb_nom(:,3))
plot3(sat_end(:,1),sat_end(:,2),sat_end(:,3))
plot3(r_sat(1:5,1),r_sat(1:5,2),r_sat(1:5,3),'ro')
legend('Nominal Satellite Orbit','Debris Orbit','Final Satellite Orbit','Burn Location')


figure; hold all
plot3(sat_nom(:,1),sat_nom(:,2),sat_nom(:,3))
plot3(sat_end(:,1),sat_end(:,2),sat_end(:,3))
legend('Nominal Satellite Orbit','Final Satellite Orbit')
