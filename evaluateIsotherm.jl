workspace();
global R_G = 8.314;
global pressure = 1e5;
global temperature = 298;
global densityAdsorbent = 1130;
global molFrac_A = 0.15;
global molFrac_B = 1-molFrac_A;

struct langmuirIsotherm
	# Component A
	qSaturationSite1_A::Float64; # Saturation capacity of site 1 [mol/kg]
	qSaturationSite2_A::Float64; # Saturation capacity of site 2 [mol/kg]
	adsorptionCoeffecientSite1_A::Float64; # Adsorption coeffecient of site 1 [m3/mol]
	adsorptionCoeffecientSite2_A::Float64; # Adsorption coeffecient of site 2 [m3/mol]
	internalEnergySite1_A::Float64; # Internal energy of site 1 [J/mol]
	internalEnergySite2_A::Float64; # Internal energy of site 2 [J/mol]

	# Component B
	qSaturationSite1_B::Float64; # Saturation capacity of site 1 [mol/kg]
	qSaturationSite2_B::Float64; # Saturation capacity of site 2 [mol/kg]
	adsorptionCoeffecientSite1_B::Float64; # Adsorption coeffecient of site 1 [m3/mol]
	adsorptionCoeffecientSite2_B::Float64; # Adsorption coeffecient of site 2 [m3/mol]
	internalEnergySite1_B::Float64; # Internal energy of site 1 [J/mol]
	internalEnergySite2_B::Float64; # Internal energy of site 2 [J/mol]
end

# Evaluate the Langmuir Isotherm either single or dual-site Langmuir
function evaluateLangmuir(materialIsotherm)
	# Calculate the concetration for species A [mol/m3]
	conc_A = pressure * molFrac_A / (R_G*temperature);

	# Site 1 loading for component A [mol/m3]
	# Adsorption coeffecient for site 1 for component A
	bsite1_A = materialIsotherm.adsorptionCoeffecientSite1_A * exp(-materialIsotherm.internalEnergySite1_A/(R_G*temperature));
	qsite1_A = (materialIsotherm.qSaturationSite1_A * densityAdsorbent * bsite1_A * conc_A) / (1 + bsite1_A * conc_A);

	# Site 2 loading for component A [mol/m3]
	# Adsorption coeffecient for site 2 for component A
	bsite2_A = materialIsotherm.adsorptionCoeffecientSite2_A * exp(-materialIsotherm.internalEnergySite2_A/(R_G*temperature));
	qsite2_A = (materialIsotherm.qSaturationSite2_A * densityAdsorbent * bsite2_A * conc_A) / (1 + bsite2_A * conc_A);

	# Equilibrium loading for component A [mol/m3]
	equilibriumLoading_A = qsite1_A + qsite2_A;

	# Calculate the concetration for species B [mol/m3]
	conc_B = pressure * molFrac_B / (R_G*temperature);

	# Site 1 loading for component B [mol/m3]
	# Adsorption coeffecient for site 1 for component B
	bsite1_B = materialIsotherm.adsorptionCoeffecientSite1_B * exp(-materialIsotherm.internalEnergySite1_B/(R_G*temperature));
	qsite1_B = (materialIsotherm.qSaturationSite1_B * densityAdsorbent * bsite1_B * conc_B) / (1 + bsite1_B * conc_B);

	# Site 2 loading for component B [mol/m3]
	# Adsorption coeffecient for site 2 for component B
	bsite2_B = materialIsotherm.adsorptionCoeffecientSite2_B * exp(-materialIsotherm.internalEnergySite2_B/(R_G*temperature));
	qsite2_B = (materialIsotherm.qSaturationSite2_B * densityAdsorbent * bsite2_B * conc_B) / (1 + bsite2_B * conc_B);

	# Equilibrium loading for component B [mol/m3]
	equilibriumLoading_B = qsite1_B + qsite2_B;

	# Return equilibrium loading for component A and B [mol/m3]
	return equilibriumLoading_A, equilibriumLoading_B
end

# Create constructor for the Langmuir Isotherm type
materialIsotherm = langmuirIsotherm(4.390,0,2.50e-06,0,-31194.15,0,
																		4.390,0,2.70e-06,0,-16375.02,0);
isothermLoading = evaluateLangmuir(materialIsotherm)./densityAdsorbent;
println(isothermLoading)
