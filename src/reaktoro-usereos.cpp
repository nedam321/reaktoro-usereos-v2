#include "reaktoro-usereos.hpp"

#include <Reaktoro/Reaktoro.hpp>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include <array>
#include <fstream>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace {

	const char* moduleName = "GEMS Test module";
	const int filenameLength = 256;
	const int titleLength = 80;

    double fluidViscosity[2] = { 1e-3, 1e-4 };
	double solidViscosity = 1e7;

} // namespace

struct Linspace {
	double from, to;
	int count;
};

std::vector<double> linspace(double begin, double end, int count) {
	std::vector<double> result(count);
	for (int i = 0; i < count; ++i) {
		double x = double(i) / (count - 1);
		result[i] = x * (end - begin) + begin;
	}
	return result;
}

std::vector<double> linspace(const Linspace& ls) { return linspace(ls.from, ls.to, ls.count); }

void to_json(json& j, const Linspace& ls) { j = json{ {"from", ls.from}, {"to", ls.to}, {"count", ls.count} }; }

void from_json(const json& j, Linspace& ls) {
	j.at("from").get_to(ls.from);
	j.at("to").get_to(ls.to);
	j.at("count").get_to(ls.count);
}

struct UserEosConfig {
	std::string databasePath;

	Linspace temperatureGrid;
	Linspace pressureGrid;

	std::string aqueousPhaseElements = "";
	std::string aqueousPhaseSpecies = "";

	std::string mineralPhaseElements = "";
	std::string mineralPhaseSpecies = "";

	std::string gaseousPhaseElements = "";
	std::string gaseousPhaseSpecies = "";
};

UserEosConfig readConfigFromFile(std::string_view path) {
    std::ifstream in(path.data());
    json j;
    in >> j;

    UserEosConfig cfg;

    j.at("database").get_to(cfg.databasePath);
    auto itp = j.at("interpolationGrid");
    itp.at("temperatureC").get_to(cfg.temperatureGrid);
    itp.at("pressureBar").get_to(cfg.pressureGrid);

    if (j.contains("aqueousPhase")) {
        auto temp = j["aqueousPhase"];
        if (temp.contains("fromElements")) {
            temp["fromElements"].get_to(cfg.aqueousPhaseElements);
        }
        if (temp.contains("fromSpecies")) {
            temp["fromSpecies"].get_to(cfg.aqueousPhaseSpecies);
        }
    }
    if (j.contains("mineralPhase")) {
        auto temp = j["mineralPhase"];
        if (temp.contains("fromElements")) {
            temp["fromElements"].get_to(cfg.mineralPhaseElements);
        }
        if (temp.contains("fromSpecies")) {
            temp["fromSpecies"].get_to(cfg.mineralPhaseSpecies);
        }
    }
    if (j.contains("gaseousPhase")) {
        auto temp = j["gaseousPhase"];
        if (temp.contains("fromElements")) {
            temp["fromElements"].get_to(cfg.gaseousPhaseElements);
        }
        if (temp.contains("fromSpecies")) {
            temp["fromSpecies"].get_to(cfg.gaseousPhaseSpecies);
        }
    }
    return cfg;
}

class UserEosModule {
public:
    UserEosModule() : m_solver(m_system), m_conds(m_system) {};//otlichie1

    int initialize(const std::string& configPath);
    void getParameters(char cmpNames[], double molWeights[], char phNames[], char auxNames[]);
    bool calculateEquilibrium(const double P, const double T, const double z[], int& nPhase, double props[],
        int8_t phaseId[], double auxArray[], int8_t mode);

    int nComponents() const;
    int nMaxPhases() const;
    int nEndMembers() const;

private:
    void updateParams();
    //void logFailedParams(const double z[], const Reaktoro::EquilibriumProblem& problem) const;
    //void logFailedParams(const double z[]) const;

    Reaktoro::ChemicalSystem m_system;
    Reaktoro::EquilibriumSolver m_solver;
    Reaktoro::EquilibriumOptions m_op;
    Reaktoro::EquilibriumConditions m_conds;

    bool a = true;

    std::vector<int> m_fluidPhasesIds;
    std::vector<int> m_solidPhasesIds;

    std::vector<std::array<char, 3>> m_componentNames;
    std::vector<std::array<char, 3>> m_phasesNames;
    std::vector<std::array<char, 8>> m_speciesNames;

    std::vector<double> m_molarWeights;
    std::vector<double> m_auxWeights;

    int m_numTransportPhases = 0;
    bool m_hasSolidPhases = false;
};

int UserEosModule::initialize(const std::string& configPath) {
    UserEosConfig cfg;

    try {
        cfg = readConfigFromFile(configPath);
    }
    catch (std::exception& e) {
        std::cerr << "Error: failed to read config file\n" << e.what() << std::endl;
    }

    std::cout << "Loaded config:\n";
    std::cout << "Database: '" << cfg.databasePath << "'\n";
    /*std::cout << "Interpolation grid:\n";*/
    /*std::cout << "  Temperature (C): [" << cfg.temperatureGrid.from << ' ' << cfg.temperatureGrid.to << "] / "
        << cfg.temperatureGrid.count << " points\n";
    std::cout << "  Pressure (bar): [" << cfg.pressureGrid.from << ' ' << cfg.pressureGrid.to << "] / "
        << cfg.pressureGrid.count << " points\n";*/

    std::cout << "Loading database... ";
    try {
        Reaktoro::AqueousPhase aq(cfg.aqueousPhaseSpecies);
        Reaktoro::GaseousPhase gas(cfg.gaseousPhaseSpecies);
        Reaktoro::MineralPhase sol(cfg.mineralPhaseSpecies);
        Reaktoro::ThermoFunDatabase db(Reaktoro::ThermoFunDatabase::fromFile(cfg.databasePath));
        m_system = Reaktoro::ChemicalSystem(db, aq, gas, sol);
    }
    catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return -1;
    }

    std::cout << "Done" << std::endl;

    m_op.warmstart = false;
    m_solver.setOptions(m_op);

    m_solver = Reaktoro::EquilibriumSolver(m_system);

    m_solver.setOptions(m_op);
    //m_conds = Reaktoro::EquilibriumConditions(m_system);
    updateParams();

    return 0;
}

void UserEosModule::updateParams() {
m_componentNames.clear();
m_phasesNames.clear();
m_speciesNames.clear();

m_molarWeights.clear();
m_auxWeights.clear();

m_fluidPhasesIds.clear();
m_solidPhasesIds.clear();

m_numTransportPhases = 0;
m_hasSolidPhases = false;

const auto& phases = m_system.phases();
for (int idx = 0; idx < phases.size(); ++idx) {
    const auto& phase = phases[idx];
    if (idx==0||idx==1) {//vot tyt drygaya proverka
    m_numTransportPhases++;
    m_fluidPhasesIds.push_back(idx);
    /*std::array<char, 3> name = {' ', ' ', ' '};
    std::copy_n(phase.name().data(), std::min(phase.name().length(), size_t(3)), name.begin());
    m_phasesNames.push_back(name);

    std::cout << "  New fluid phase: " << phase.name() << std::endl;*/
    if (idx == 0) {
        std::cout << "  New fluid phase: LIQ" << std::endl;
        m_phasesNames.push_back({ 'L', 'I', 'Q' });
    }
    if (idx == 1) {
        std::cout << "  New fluid phase: GAS" << std::endl;                                                              
        m_phasesNames.push_back({ 'G', 'A', 'S' });
    }
    }
    else if (idx==2) {//tyt tozhe pomenyat proverky SoM ili agrstate
         if (!m_hasSolidPhases) {
             m_numTransportPhases++;
             m_hasSolidPhases = true;
         }
         m_solidPhasesIds.push_back(idx);
     }
}

if (m_hasSolidPhases) {
     std::cout << "  New solid phase: SOL" << std::endl;
     m_phasesNames.push_back({ 'S', 'O', 'L' });
 }
std::array<char, 3> name = { '1', ' ', ' ' };
m_componentNames.push_back(name);
name = { '2', ' ', ' ' };
m_componentNames.push_back(name);
name = { '3', ' ', ' ' };
m_componentNames.push_back(name);
name = { '4', ' ', ' ' };
m_componentNames.push_back(name);
m_auxWeights.push_back(0.0180153);
m_auxWeights.push_back(0.0440096);
m_auxWeights.push_back(0.100087);
m_auxWeights.push_back(0.100087);
m_molarWeights.push_back(0.0180153);
m_molarWeights.push_back(0.0440096);
m_molarWeights.push_back(0.100087);
m_molarWeights.push_back(0.100087);
std::cout << "Elements: ";
/*for (const auto& elem : m_system.elements()) {
    //std::array<char, 3> name = { '1', ' ', ' ' };
    //std::copy_n(elem.name().data(), std::min(elem.name().length(), size_t(3)), name.begin());
    //m_componentNames.push_back(name);
    m_molarWeights.push_back(elem.molarMass());
    std::cout << elem.name() << ' ';
}*/
std::cout << std::endl;

std::cout << "Species: ";
for (const auto& specie : m_system.species()) {
    std::array<char, 8> name = { ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ' };
    std::string aux("#");
    std::string esc = specie.name().substr(0, 7);
    aux.append(esc);
    //std::cout << specie.molarMass() << " " << aux << std::endl;
    std::copy_n(aux.data(), std::min(aux.length(), size_t(8)), name.begin());
    m_speciesNames.push_back(name);

    std::cout << specie.name() << ' ';
}
std::cout << std::endl;
}

void UserEosModule::getParameters(char cmpNames[], double molWeights[], char phNames[], char auxNames[]) {
    memcpy(cmpNames, m_componentNames.data(), 3 * size_t(nComponents()));
    memcpy(phNames, m_phasesNames.data(), 3 * size_t(nMaxPhases()));
    if (auxNames != nullptr) {
        memcpy(auxNames, m_speciesNames.data(), 8 * size_t(nEndMembers()));
    }
    std::copy(m_molarWeights.begin(), m_molarWeights.end(), molWeights);
}

bool UserEosModule::calculateEquilibrium(const double P, const double T, const double z[], int& nPhase, double props[],
    int8_t phaseId[], double auxArray[], int8_t mode) {
    Reaktoro::ChemicalState state(m_system);
    Reaktoro::EquilibriumConditions m_conds(m_system);
    //std::cout << "P=" << P << " T=" << T << std::endl;
    //for (int ic = 0; ic < nComponents(); ++ic)
        //std::cout << z[ic] << " ";
    //std::cout << std::endl;

    m_conds.pressure(P, "Pa");
    m_conds.temperature(T, "K");

    Reaktoro::VectorXd elems = Reaktoro::zeros(m_system.elements().size() + 1);
    elems[0] = z[0] / m_auxWeights[0] * 2.;
    //std::cout << elems[0]<<" ";
    elems[1] = z[1] / m_auxWeights[1] + z[2] / m_auxWeights[2] + z[3] / m_auxWeights[3] * 0.5;
    //std::cout << elems[1] << " ";
    elems[2] = z[0] / m_auxWeights[0] + z[1] / m_auxWeights[1] * 2. + z[2] / m_auxWeights[2] * 3. + z[3] / m_auxWeights[3] * 1.5;
    //std::cout << elems[2] << " ";
    elems[3] = z[2] / m_auxWeights[2] + z[3] / m_auxWeights[3];
    //std::cout << elems[3] << " ";
    //std::cout << std::endl;
    /*double summ = 0;
    for (int i = 0; i < 4; i++)
        summ += elems[i];
    std::cout << summ << std::endl;*/
    m_conds.setInitialComponentAmounts(elems);

    Reaktoro::EquilibriumResult result;

    /* try {
         m_solver.approximate(state, m_problem);
     }
     catch (std::exception& e) {
         std::cerr << "Error (approximate): " << e.what() << std::endl;
         std::cerr << "System params:" << std::endl;
         logFailedParams(z);
         state.output("failed.txt");
         exit(1);
     }*/

    m_solver = Reaktoro::EquilibriumSolver(m_system);

    try {
        result = m_solver.solve(state, m_conds);
            //state.output(std::cout);
            //std::cout << std::endl;
    }
    catch (std::exception& e) {
        std::cerr << "Error (solve): " << e.what() << std::endl;
        std::cerr << "System params:" << std::endl;
        //logFailedParams(z);
        state.output("failed.txt");
        exit(1);
    }

    // if (!result.optimum.succeeded) {
    //   std::cerr << "Error (convergence check): the problem did not converge, system params:" << std::endl;
    //   logFailedParams(z);
    //   state.output("failed.txt");
    //   exit(1);
    // }
    //state.output(std::cout);
    //std::cout << std::endl;
    size_t offset = size_t(6 + nComponents());
    std::fill(props, props + nMaxPhases() * offset, 0.0);

    Reaktoro::ChemicalProps properties = state.props();
    auto aqprops = properties.phaseProps(m_fluidPhasesIds[0]);
    auto phaseMasses = aqprops.mass();
    auto phaseVolumes = aqprops.volume();
    auto phaseDensities = aqprops.density();
    auto phaseEnthalpies = aqprops.specificEnthalpy();
    //auto phaseMasses = properties.phaseMasses();
    //auto phaseVolumes = properties.phaseVolumes();
    //auto phaseDensities = properties.phaseDensities();
    //auto phaseEnthalpies = properties.phaseSpecificEnthalpies();

    auto totalVolume = double(properties.volume());

    for (int ip = 0; ip < m_fluidPhasesIds.size(); ip++) {
        int iphase = m_fluidPhasesIds[ip];
        auto aqprops = properties.phaseProps(iphase);
        auto phaseMasses = aqprops.mass();
        auto phaseVolumes = aqprops.volume();
        auto phaseDensities = aqprops.density();
        auto phaseEnthalpies = aqprops.specificEnthalpy();

        //props[ip * offset + 0] = double(phaseDensities[iphase]);
        //props[ip * offset + 1] = double(phaseEnthalpies[iphase]);
        //props[ip * offset + 2] = fluidViscosity;
        //props[ip * offset + 3] = double(phaseVolumes[iphase]) / totalVolume;
        //props[ip * offset + 4] = double(phaseVolumes[iphase]) / totalVolume;

        props[ip * offset + 0] = double(phaseDensities);
        props[ip * offset + 1] = double(phaseEnthalpies);
        props[ip * offset + 2] = fluidViscosity[ip];
        props[ip * offset + 3] = double(phaseVolumes) / totalVolume;
        props[ip * offset + 4] = double(phaseVolumes) / totalVolume;
        //if (ip == 0)
            //std::cout << " LIQ:     ";
        //if (ip == 1)
            //std::cout << " GAS:     ";
        //for (int i = 0; i < 5;i++)
            //std::cout << props[ip * offset + i] << " ";
        //std::cout << std::endl;

        double sumC = 0.0;
        if (ip == 0) {
            double x1 = properties.elementAmountInPhase(0, iphase);
            props[ip * offset + 6 + 0] = x1 / 2 * m_auxWeights[0] / phaseMasses;
            double x2 = properties.elementAmountInPhase(3, iphase);
            props[ip * offset + 6 + 3] = x2 * m_auxWeights[3] / phaseMasses;
            double x3 = properties.elementAmountInPhase(1, iphase);
            if(x3<=x2)
                props[ip * offset + 6 + 1] = 0;
            else
                props[ip * offset + 6 + 1] = (x3 - x2) * m_auxWeights[1] / phaseMasses;
            props[ip * offset + 6 + 2] = 0;
        }
        if (ip == 1) {
            double x1 = properties.elementAmountInPhase(0, iphase);
            props[ip * offset + 6 + 0] = x1 / 2 * m_auxWeights[0] / phaseMasses;
            double x2 = properties.elementAmountInPhase(1, iphase);
            props[ip * offset + 6 + 1] = x2 * m_auxWeights[1] / phaseMasses;
                props[ip * offset + 6 + 2] = 0;
            props[ip * offset + 6 + 3] = 0;
        }
        for (int ic = 0; ic < nComponents(); ++ic) {
            //double x = properties.elementAmountInPhase(ic, iphase);
            //props[ip * offset + 6 + ic] = x * m_molarWeights[ic] / phaseMasses;
            //std::cout << props[ip * offset + 6 + ic] << " ";
            sumC += props[ip * offset + 6 + ic];
            std::cout << props[ip * offset + 6 + ic] << "  ";
        }
        std::cout << std::endl;

        if (std::abs(sumC - 1.0) > 1e-12) {
            std::cerr << "Error: fluid concentrations do not sum to 1!!! "<<sumC<<"\n";
        }

        phaseId[ip] = int8_t(ip + 1);
    }

    nPhase = nMaxPhases();

    if (m_hasSolidPhases) {
        double solidVolume = 0.0;
        int ip = m_fluidPhasesIds.size();

        double solidMass = 0.0;
        double solidEnthalpy = 0.0;

        for (int is = 0; is < m_solidPhasesIds.size(); ++is) {
            int iphase = m_solidPhasesIds[is];
            auto aqprops = properties.phaseProps(iphase);
            auto phaseMasses = aqprops.mass();
            auto phaseVolumes = aqprops.volume();
            auto phaseDensities = aqprops.density();
            auto phaseEnthalpies = aqprops.specificEnthalpy();
            solidMass += double(phaseMasses);
            solidEnthalpy += double(phaseEnthalpies) * double(phaseMasses);
            solidVolume += phaseVolumes;

            //for (int ic = 0; ic < nComponents(); ++ic) {
                double x = properties.elementAmountInPhase(3, iphase);
                props[ip * offset + 6 + 2] += x * m_auxWeights[3];
            //}
        }

        if (solidMass < 1e-12) {
            nPhase--;

            props[ip * offset + 0] = NAN;
            props[ip * offset + 1] = NAN;
            props[ip * offset + 2] = NAN;
            props[ip * offset + 3] = NAN;
            props[ip * offset + 4] = NAN;
            //std::cout << "SOL          ";
            //for (int i = 0; i < 5; i++)
                //std::cout << props[ip * offset + i] << " ";
            //std::cout << std::endl;
        }
        else {
            // Specific enthalpy
            solidEnthalpy /= solidMass;

            props[ip * offset + 0] = solidMass / solidVolume;
            props[ip * offset + 1] = solidEnthalpy;
            props[ip * offset + 2] = solidViscosity;
            props[ip * offset + 3] = solidVolume / totalVolume;
            props[ip * offset + 4] = 0.0;
            //std::cout << "SOL          ";
            //for(int i=0;i<5;i++)
            //std::cout << props[ip * offset + i] << " ";
            //std::cout << std::endl;
            double sumC = 0.0;
            // Normalize the concentrations
            for (int ic = 0; ic < nComponents(); ++ic) {
                props[ip * offset + 6 + ic] /= solidMass;
                //std::cout << props[ip * offset + 6 + ic] << " ";
                sumC += props[ip * offset + 6 + ic];
            }
            //std::cout << std::endl;
           if (std::abs(sumC - 1.0) > 1e-12) {
                std::cerr << "Error: solid concentrations do not sum to 1!!! "<<sumC<<"\n";
            }
            phaseId[ip] = int8_t(ip + 1);
        }
    }

    for (int iem = 0; iem < nEndMembers(); ++iem) {
        auxArray[iem] = state.speciesAmount(iem); // моли возвращает
    }

    return true;
}

int UserEosModule::nComponents() const { return m_system.elements().size(); }

int UserEosModule::nMaxPhases() const { return m_numTransportPhases; }

int UserEosModule::nEndMembers() const { return m_system.species().size(); }

static UserEosModule userModule;

void ReadConfigurationFile(char* filename, int* ierr, char* title) {
    std::string cFilename(filename, filename + filenameLength);
    std::cout << cFilename << std::endl;
    size_t spacePos = cFilename.find(' ');
    if (spacePos != std::string::npos) {
        cFilename[cFilename.find(' ')] = '\0';
    }
    memset(title, ' ', size_t(titleLength));
    memcpy(title, moduleName, strlen(moduleName));

    try {
        *ierr = userModule.initialize(cFilename);
    }
    catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

void GetDimensions(int* nComponents, int* nPhaseMax, int* nAux) {
    *nComponents = userModule.nComponents();
    *nPhaseMax = userModule.nMaxPhases();
    *nAux = userModule.nEndMembers();

    std::cout << "nMaxPhases: " << userModule.nMaxPhases() << std::endl;
}

void GetGlobalParameters(char* cmpNames, double* molWeights, char* phNames, char* auxNames, char* auxUnits,
    int8_t* opt) {
    try {
        userModule.getParameters(cmpNames, molWeights, phNames, auxNames);
    }
    catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    char unit[] = "NODIM   ";
    for (int idx = 0; idx < userModule.nEndMembers(); ++idx) {
        memcpy(auxUnits + size_t(idx) * 8, unit, 8);
    }
    opt[0] = 1;
    opt[1] = 0;

}

void PhaseEquilibrium(double* pres, double* temp, double* z, int* nPhase, double* props, int8_t* phaseId,
    double* auxArray, int8_t* mode) {
    userModule.calculateEquilibrium(*pres, *temp, z, *nPhase, props, phaseId, auxArray, *mode);
}
