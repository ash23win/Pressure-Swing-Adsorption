#workspace();
using DifferentialEquations

tic();

flagPlotIsotherm = true; # Flag for plotting the isotherms

# Sturcture for isotherm properties
struct langmuirIsotherm
	# Component A
	qSaturationSite1_A::Float64; # Saturation capacity of site 1 [mol/m3]
	qSaturationSite2_A::Float64; # Saturation capacity of site 2 [mol/m3]
	adsorptionCoeffecientSite1_A::Float64; # Adsorption coeffecient of site 1 [m3/mol]
	adsorptionCoeffecientSite2_A::Float64; # Adsorption coeffecient of site 2 [m3/mol]
	internalEnergySite1_A::Float64; # Internal energy of site 1 [J/mol]
	internalEnergySite2_A::Float64; # Internal energy of site 2 [J/mol]

	# Component B
	qSaturationSite1_B::Float64; # Saturation capacity of site 1 [mol/m3]
	qSaturationSite2_B::Float64; # Saturation capacity of site 2 [mol/m3]
	adsorptionCoeffecientSite1_B::Float64; # Adsorption coeffecient of site 1 [m3/mol]
	adsorptionCoeffecientSite2_B::Float64; # Adsorption coeffecient of site 2 [m3/mol]
	internalEnergySite1_B::Float64; # Internal energy of site 1 [J/mol]
	internalEnergySite2_B::Float64; # Internal energy of site 2 [J/mol]
end

# CONSTANTS
const UNIVERSAL_GAS_CONSTANT = 8.314472; # Universal gas constant [------]
const NUMBER_OF_COMPONENTS = 2; # Number of components in the feed gas
const NUMBER_OF_GRID_POINTS = 30; # Number of grid points along spatial discretization
const MEMORY_ALLOCATION_POINTS = 10000; # Number of points to be allocated for the main array

# FEED/PROCESS CONDITIONS
const PROCESS_FEEDTEMPERATURE = 298.15; # Temperature of the feed [K]
const PROCESS_AMBIENTTEMPERATURE = 298.15; # Temperature of the ambient atmosphere [K]
const PROCESS_HIGHPRESSURE = 1e5; # High pressure of the process [Pa]
const PROCESS_INTPRESSURE = 0.10e5; # Intermediate pressure of the process [Pa]
const PROCESS_LOWPRESSURE = 0.03e5; # Low pressure of the process [Pa]
const PROCESS_TIMEASDORPTION = 50; # Time of the adsorption step [s]
const PROCESS_FEEDVELOCITY = 1; # Feed velocity [m/s]
const PROCESS_FEEDCONC_A = 0.15; # Feed concentration of component A [-]
const PROCESS_FEEDCONC_B = 1-PROCESS_FEEDCONC_A; # Feed concentration of component B [-]
const PROCESS_INITIALCONC_A = 0; # Initial gas phase concentration of component A [-]
const PROCESS_INITIALCONC_B = 1-PROCESS_INITIALCONC_A; # Initial gas phase concentration of component A [-]

# REFERENCE CONDITIONS
const REFERENCE_PRESSURE = PROCESS_HIGHPRESSURE; # Reference pressure [Pa]
const REFERENCE_TEMPERATURE = PROCESS_FEEDTEMPERATURE; # Reference temperature [K]

# MIXTURE PROPERTIES
const GAS_SPECIFICHEATCAPACITY = 1010.6; # Specific heat capacity of the gas mixture [J/kg/K]
const GAS_massTransferCoeff_A = 0.1; # Mass transfer coeffecient of component A [-]
const GAS_massTransferCoeff_B = 0.1; # Mass transfer coeffecient of component B [-]
const GAS_MOLECULARDIFFUSIVITY = 1.2955e-5; #Molecular diffusivity of CO2-N2 mixture at 298.15K [m2/s]

# ADSORBENT PROPERTIES
const ADSORBENT_DENSITY = 1130; # Density of the adsorbent [kg/m3]
const ADSORBENT_PARTICLERADIUS = 0.00075; # Adsorbent particle radius [m]
const ADSORBENT_SPECIFICHEATCAPACITY = 1070.0; # Specific heat capacity of the adsorbent [J/kg/K]

# COLUMN PROPERTIES
const COLUMN_DENSITY = 7800; # Density of the column wall [kg/m2]
const COLUMN_BEDLENGTH = 1.00; # Length of the column [m]
const COLUMN_INNERRADIUS = 0.1445 # Inner radius of the column [m]
const COLUMN_OUTERRADIUS = 0.162 # Inner radius of the column [m]
const COLUMN_BEDVOIDAGE = 0.37; # Porosity of the bed [-]
const COLUMN_TORTUOSITY = 3.0; # Tortuosity [-]
const COLUMN_TEMPERATURE = 298.15; # Temperature of the ambient atmosphere [K]
const COLUMN_THERMALCONDUCTIVITY = 16.0; # Thermal conducitivty of the column wall [W/m/K]
const COLUMN_INNERHEATCOEFFECIENT = 8.6; # Heat transfer coefficient inside the column [W/m2/K]
const COLUMN_OUTERHEATCOEFFECIENT = 2.5; # Heat transfer coefficient outside the column [W/m2/K]
const COLUMN_HEATCONDUCTIONCOEFFECIENT = 0.0903; # Effective heat conduction coeffecient [W/m/K]
const COLUMN_SPECIFICHEATCAPACITY = 502; # Specific heat capacity of the column wall [J/kg/K]
const COLUMN_AXIALDISPERSION = 0.7*(GAS_MOLECULARDIFFUSIVITY*REFERENCE_PRESSURE/PROCESS_HIGHPRESSURE) + 0.5*PROCESS_FEEDVELOCITY*(2.0*ADSORBENT_PARTICLERADIUS);
const COLUMN_MACROPOREDIFFUSIVITY = GAS_MOLECULARDIFFUSIVITY/COLUMN_TORTUOSITY; # Macropore diffusivity [m2/s]

# NON-DIMENSIONAL QUANTITIES
const NONDIMENSIONAL_TIME = COLUMN_BEDLENGTH/PROCESS_FEEDVELOCITY; # Non-dimensional time [-]
const NONDIMENSIONAL_DELZ = 1/NUMBER_OF_GRID_POINTS; # Length of the cell along spatial discretization [-]
const NONDIMENSIONAL_PECLETMASS = PROCESS_FEEDVELOCITY*COLUMN_BEDLENGTH/COLUMN_AXIALDISPERSION; # Axial dispersion along the column

# ADSORBENT PROPERTIES
# Component A
qSat1_A = 4.390*ADSORBENT_DENSITY;
qSat2_A = 0.000*ADSORBENT_DENSITY;

site1b0_A = 2.50e-06;
site2b0_A = 0;

site1U_A = -31194.15;
site2U_A = 0;

# Component B
qSat1_B = 4.390*ADSORBENT_DENSITY;
qSat2_B = 0.000*ADSORBENT_DENSITY;

site1b0_B = 2.70e-06;
site2b0_B = 0;

site1U_B = -16375.02;
site2U_B = 0;

# Reference saturation capacity
qSat_REF = qSat1_B;

# Heats of adsorption for component A and B for individual sites [J/mol]
site1H_A = site1U_A - UNIVERSAL_GAS_CONSTANT*REFERENCE_TEMPERATURE;
site2H_A = site2U_A - UNIVERSAL_GAS_CONSTANT*REFERENCE_TEMPERATURE;
site1H_B = site1U_B - UNIVERSAL_GAS_CONSTANT*REFERENCE_TEMPERATURE;
site2H_B = site2U_B - UNIVERSAL_GAS_CONSTANT*REFERENCE_TEMPERATURE;

# Create structure using single site Langmuir Isotherm
materialIsotherm = langmuirIsotherm(qSat1_A/qSat_REF,qSat2_A/qSat_REF,site1b0_A,site2b0_A,site1U_A,site2U_A,
																		qSat1_B/qSat_REF,qSat2_B/qSat_REF,site1b0_B,site2b0_B,site1U_B,site2U_B);

# Heats of adsorption for component A and B obtained by combining the sites using a reference
delH_A = qSat1_A/qSat_REF*site1H_A + qSat2_A/qSat_REF*site2H_A;
delH_B = qSat1_B/qSat_REF*site1H_B + qSat2_B/qSat_REF*site2H_B;

# Vector of heats of adsorption
delH = [delH_A, delH_B];

# Evaluate the Langmuir Isotherm either single or dual-site Langmuir
function evaluateLangmuir(materialIsotherm,currentPressure,currentTemperature,currentMoleFraction_A)
	# Calculate the concetration for species A [mol/m3]
	conc_A = currentPressure * currentMoleFraction_A / (UNIVERSAL_GAS_CONSTANT*currentTemperature);

	# Site 1 loading for component A [mol/m3]
	# Adsorption coeffecient for site 1 for component A
	bsite1_A = materialIsotherm.adsorptionCoeffecientSite1_A * exp(-materialIsotherm.internalEnergySite1_A/(UNIVERSAL_GAS_CONSTANT*currentTemperature));
	qsite1_A = (materialIsotherm.qSaturationSite1_A * bsite1_A * conc_A) / (1 + bsite1_A * conc_A);

	# Site 2 loading for component A [mol/m3]
	# Adsorption coeffecient for site 2 for component A
	bsite2_A = materialIsotherm.adsorptionCoeffecientSite2_A * exp(-materialIsotherm.internalEnergySite2_A/(UNIVERSAL_GAS_CONSTANT*currentTemperature));
	qsite2_A = (materialIsotherm.qSaturationSite2_A * bsite2_A * conc_A) / (1 + bsite2_A * conc_A);

	# Equilibrium loading for component A [mol/m3]
	equilibriumLoading_A = qsite1_A + qsite2_A;

	# Calculate the concetration for species B [mol/m3]
	conc_B = currentPressure * (1 - currentMoleFraction_A) / (UNIVERSAL_GAS_CONSTANT*currentTemperature);

	# Site 1 loading for component B [mol/m3]
	# Adsorption coeffecient for site 1 for component B
	bsite1_B = materialIsotherm.adsorptionCoeffecientSite1_B * exp(-materialIsotherm.internalEnergySite1_B/(UNIVERSAL_GAS_CONSTANT*currentTemperature));
	qsite1_B = (materialIsotherm.qSaturationSite1_B * bsite1_B * conc_B) / (1 + bsite1_B * conc_B);

	# Site 2 loading for component B [mol/m3]
	# Adsorption coeffecient for site 2 for component B
	bsite2_B = materialIsotherm.adsorptionCoeffecientSite2_B * exp(-materialIsotherm.internalEnergySite2_B/(UNIVERSAL_GAS_CONSTANT*currentTemperature));
	qsite2_B = (materialIsotherm.qSaturationSite2_B * bsite2_B * conc_B) / (1 + bsite2_B * conc_B);

	# Equilibrium loading for component B [mol/m3]
	equilibriumLoading_B = qsite1_B + qsite2_B;

	# Return equilibrium loading for component A and B [mol/m3]
	return equilibriumLoading_A, equilibriumLoading_B
end

function runAdsorption(t,mainArrayTemp,mainArrayDerivative)
	nCount = 0;
	velocityAdsorption = 1.0;
	NONDIMENSIONAL_DELT = 1;
	# Evaluate the equilibrium using the isotherm
	isothermLoadings = evaluateLangmuir.(materialIsotherm,mainArrayTemp[5*NUMBER_OF_GRID_POINTS+1:6*NUMBER_OF_GRID_POINTS]*REFERENCE_PRESSURE,mainArrayTemp[4*NUMBER_OF_GRID_POINTS+1:5*NUMBER_OF_GRID_POINTS]*REFERENCE_TEMPERATURE,mainArrayTemp[1:NUMBER_OF_GRID_POINTS]);

	isothermLoading_A = zeros(NUMBER_OF_GRID_POINTS);
	isothermLoading_B = zeros(NUMBER_OF_GRID_POINTS);

	for i=1:length(isothermLoading_A)
		isothermLoading_A[i] = isothermLoadings[i][1];
		isothermLoading_B[i] = isothermLoadings[i][2];
	end

	# Evaluate the Linear Driving Force model for Component A
	mainArrayDerivative[1*NUMBER_OF_GRID_POINTS+1:2*NUMBER_OF_GRID_POINTS] = GAS_massTransferCoeff_A*(isothermLoading_A - mainArrayTemp[NUMBER_OF_GRID_POINTS+1:2*NUMBER_OF_GRID_POINTS]);

	# Evaluate the Linear Driving Force model for Component B
	mainArrayDerivative[2*NUMBER_OF_GRID_POINTS+1:3*NUMBER_OF_GRID_POINTS] = GAS_massTransferCoeff_B*(isothermLoading_B - mainArrayTemp[2*NUMBER_OF_GRID_POINTS+1:3*NUMBER_OF_GRID_POINTS]);
  #
	# moleFraction_A_IN = (PROCESS_FEEDCONC_A*velocityAdsorption*NONDIMENSIONAL_PECLETMASS*NONDIMENSIONAL_DELZ/2 + mainArray[1,nCount])/(velocityAdsorption*NONDIMENSIONAL_PECLETMASS*NONDIMENSIONAL_DELZ/2);
end

###### MAIN FUNCTION ######
## INITIAL CONDITIONS FOR THE BED
# Initialize main array as an empty array
mainArray = zeros(6*NUMBER_OF_GRID_POINTS,1);

# COMPONENT A MOLE FRACTION IN GAS PHASE [1:NUMBER_OF_GRID_POINTS]
mainArray[0*NUMBER_OF_GRID_POINTS+1:1*NUMBER_OF_GRID_POINTS] = ones(1,NUMBER_OF_GRID_POINTS).*PROCESS_FEEDCONC_A;

# INITIALIZE THE SOLID PHASE AT INITIAL CONDITIONS CONDITIONS
# Evaluate the Langmuir isotherm for the two components
isothermLoadings = evaluateLangmuir(materialIsotherm,PROCESS_HIGHPRESSURE,PROCESS_FEEDTEMPERATURE,PROCESS_INITIALCONC_A);

# Component A
mainArray[1*NUMBER_OF_GRID_POINTS+1:2*NUMBER_OF_GRID_POINTS] = ones(1,NUMBER_OF_GRID_POINTS).*isothermLoadings[1][1];
# Component B
mainArray[2*NUMBER_OF_GRID_POINTS+1:3*NUMBER_OF_GRID_POINTS] = ones(1,NUMBER_OF_GRID_POINTS).*isothermLoadings[2][1];

# INITIALIZE THE BED TEMPERATURE WITH THE FEED TEMPERATURE
mainArray[3*NUMBER_OF_GRID_POINTS+1:4*NUMBER_OF_GRID_POINTS] = ones(1,NUMBER_OF_GRID_POINTS).*PROCESS_FEEDTEMPERATURE/REFERENCE_TEMPERATURE;

# INITIALIZE THE COLUMN TEMPERATURE WITH THE FEED TEMPERATURE
mainArray[4*NUMBER_OF_GRID_POINTS+1:5*NUMBER_OF_GRID_POINTS] = ones(1,NUMBER_OF_GRID_POINTS).*COLUMN_TEMPERATURE/REFERENCE_TEMPERATURE;

# INITIALIZE THE COLUMN PRESSURE WITH THE FEED PRESSURE
mainArray[5*NUMBER_OF_GRID_POINTS+1:6*NUMBER_OF_GRID_POINTS,1] = ones(1,NUMBER_OF_GRID_POINTS).*PROCESS_HIGHPRESSURE/REFERENCE_PRESSURE;

# SOLVE THE ODE
tspan = (0.0,20.0);
prob = ODEProblem(runAdsorption,mainArray,tspan)
sol = solve(prob)

if flagPlotIsotherm
	pressureValues = 0:1e4:1e5;
	tempValues = ones(size(pressureValues)).*298.15;
	molValues = ones(size(pressureValues)).*0.15;
	isothermLoadings = evaluateLangmuir.(materialIsotherm,pressureValues,tempValues,molValues);
	isothermLoading_A = zeros(size(pressureValues));
	isothermLoading_B = zeros(size(pressureValues));
	for i=1:length(isothermLoading_A)
		isothermLoading_A[i] = isothermLoadings[i][1];
		isothermLoading_B[i] = isothermLoadings[i][2];
	end
	# gui()
	# using PyPlot
	# plot(pressureValues,isothermLoading_A)
	# plot(pressureValues,isothermLoading_B)
end
toc()
