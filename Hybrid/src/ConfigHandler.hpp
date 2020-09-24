#ifndef CONFIG_HANDLER
#define CONFIG_HANDLER

#include <iomanip>
#include <cstdlib>
#include <iostream>
//#include <libconfig.h++>
#include <string>
#include "Util.hpp"

//using namespace libconfig;

class ConfigHandler
{
private:
    //Config cfg;
    bool created_ = true;

    void getInputValues()
    { 
        /* const Setting& conf = this->cfg.getRoot()["agents"]["input"];

        try
        { */
            this->input.seed = 7;
            this->input.timeMax = 96;
            this->input.initialCondition = "input/ABM_INPUT.dat";
            this->input.fileFormat = BI_D; //2D
            //this->input.fileFormat = TRI_D; //3D
        /* }
        catch(const SettingNotFoundException &nfex)
        {
            this->created_ = false;
            std::cerr << "Input Settings Not Found!" << std::endl;
        } */
    }

    void getOutputValues()
    {
        /* const Setting& conf = this->cfg.getRoot()["agents"]["output"];

        try
        { */
            this->output.paths.agent = "output/";
            this->output.paths.nut = "output/nut/";
            this->output.filenames.number = 1;
            this->output.filenames.agent = "out";
            this->output.filenames.nut = "nut";
            this->output.nut = false;
            this->output.files = false;
            this->output.prints = false;

            this->output.justLastFile = false;
            this->output.justLastPrint = false;
        /* }
        catch(const SettingNotFoundException &nfex)
        {
            this->created_ = false;
            std::cerr << "Output Settings Not Found!" << std::endl;
        } */
    }

    void getContinuumValues()
    {
        /* const Setting& conf = this->cfg.getRoot()["agents"]["continuum"];
        try
        { */
            this->continuum.oxgD = 4.64894;
            this->continuum.oConsumptionBg = 0.0;
            this->continuum.oBorder = 1.0;

            this->continuum.hCoarse = 10.0;
            this->continuum.hRefined = 1.0;
            this->continuum.deltaT = 1.0;
        /* }
        catch(const SettingNotFoundException &nfex)
        {
            this->created_ = false;
            std::cerr << "Continuum Settings Not Found!" << std::endl;
        } */
    }

    void getAgentValues()
    {
        /* const Setting& conf = this->cfg.getRoot()["agents"]["agent"];
        try
        { */
            this->agent.nucleusRadius = 5.295;
            this->agent.radius = 9.953;
            this->agent.actionRadius = 1.214;
            this->agent.oConsumption = 0.06588;
            this->agent.actionRadius *= this->agent.radius;
        /* }
        catch(const SettingNotFoundException &nfex)
        {
            this->created_ = false;
            std::cerr << "Agent Settings Not Found!" << std::endl;
        } */
    }

    void getForcesValues()
    {
        /* const Setting& conf = this->cfg.getRoot()["agents"]["forces"];
        try
        { */
            this->forces.c_cca = -0.488836;
            this->forces.c_ccr = -10.0;
            this->forces.K = 0.1;
            this->forces.c_ct = -10.0;//*K
            this->forces.c_rct = -4.88836;//*K
            this->forces.c_ct *= this->forces.K;
            this->forces.c_rct *= this->forces.K;

        /* }
        catch(const SettingNotFoundException &nfex)
        {
            this->created_ = false;
            std::cerr << "Forces Settings Not Found!" << std::endl;
        } */
    }

    void getParametersValues()
    {
        /* const Setting& conf = this->cfg.getRoot()["agents"]["parameters"];
        try
        { */
            this->parameters.tauA = 8.6;
            this->parameters.tauP = 18.0;
            this->parameters.tauG1 = 9.0;
            this->parameters.delta_tt = 1.0;
            this->parameters.sigmaH = 0.14681;
            this->parameters.gammaA = 100.0;
            this->parameters.alphaA = 0.0;
            this->parameters.alphaP = 1.0;
        /* }
        catch(const SettingNotFoundException &nfex)
        {
            this->created_ = false;
            std::cerr << "Parameters Settings Not Found!" << std::endl;
        } */
    }

public:
    ConfigHandler (std::string configFile = "config.cfg")
    {
        /* this->created_ = false;
        try
        {
            this->cfg.readFile(configFile.c_str());
            this->created_ = true;
        }
        catch(const FileIOException &fioex)
        {
            std::cerr << "I/O error while reading file." << std::endl;
        }
        catch(const ParseException &pex)
        {
            std::cerr   << "Parse error at " << pex.getFile() << ":" << pex.getLine()
                        << " - " << pex.getError() << std::endl;
        } */
        
        this->getInputValues();
        this->getOutputValues();
        this->getAgentValues();
        this->getContinuumValues();
        this->getForcesValues();
        this->getParametersValues();
    }

    bool created(){
        return this->created_;
    }

    struct
    {
        std::string initialCondition;
        int seed;
        int timeMax;
        FileFormat fileFormat;
    }input;

    struct
    {
        struct
        {
            std::string agent;
            std::string nut;
        } paths;
        struct
        {
            std::string agent;
            std::string nut;
            int number;
        }filenames;
        bool nut;
        bool files;
        bool prints;
        bool justLastFile;
        bool justLastPrint;
    } output;

    struct
    {
        double oxgD;
        double oConsumptionBg;
        double oBorder;
        double hCoarse;
        double hRefined;
        double deltaT;
    }continuum;

    struct
    {
        double nucleusRadius;
        double radius;
        double actionRadius;
        double oConsumption;
    }agent;

    struct
    {
        double c_cca;
        double c_ccr;
        double K;
        double c_ct;
        double c_rct;
    }forces;

    struct
    {
        double tauA;
        double tauP;
        double tauG1;
        double delta_tt;
        double alphaP;
        double sigmaH;
        double alphaA;
        double gammaA;
    }parameters;

};

#endif /* end of include guard: CONFIG_HANDLER */
