#ifndef CONFIG_HANDLER
#define CONFIG_HANDLER

#include <iomanip>
#include <cstdlib>
#include <iostream>
#include <libconfig.h++>
#include <string>
#include "Util.hpp"

using namespace libconfig;

class ConfigHandler
{
private:
    Config cfg;
    bool created_;

    void getInputValues()
    {
        const Setting& conf = this->cfg.getRoot()["agents"]["input"];

        try
        {
            conf.lookupValue("seed", this->input.seed);
            conf.lookupValue("time-max", this->input.timeMax);
            conf.lookupValue("initial-condition", this->input.initialCondition);
            std::string fileFormat;
            conf.lookupValue("file-format", fileFormat);
            if (fileFormat == "2d" || fileFormat == "2D") {
                this->input.fileFormat = BI_D;
            } else if (fileFormat == "3d" || fileFormat == "3D") {
                this->input.fileFormat = TRI_D;
            }
        }
        catch(const SettingNotFoundException &nfex)
        {
            this->created_ = false;
            std::cerr << "Input Settings Not Found!" << std::endl;
        }
    }

    void getOutputValues()
    {
        const Setting& conf = this->cfg.getRoot()["agents"]["output"];

        try
        {
            conf["paths"].lookupValue("agent", this->output.paths.agent);
            conf["paths"].lookupValue("nut", this->output.paths.nut);
            conf["filenames"].lookupValue("number", this->output.filenames.number);
            conf["filenames"].lookupValue("agent", this->output.filenames.agent);
            conf["filenames"].lookupValue("nut", this->output.filenames.nut);
            conf.lookupValue("nut", this->output.nut);
            conf.lookupValue("files", this->output.files);
            conf.lookupValue("prints", this->output.prints);

            conf.lookupValue("just-last-file", this->output.justLastFile);
            conf.lookupValue("just-last-print", this->output.justLastPrint);
        }
        catch(const SettingNotFoundException &nfex)
        {
            this->created_ = false;
            std::cerr << "Output Settings Not Found!" << std::endl;
        }
    }

    void getContinuumValues()
    {
        const Setting& conf = this->cfg.getRoot()["agents"]["continuum"];
        try
        {
            conf.lookupValue("oxg-D", this->continuum.oxgD);
            conf.lookupValue("o-consumption-bg", this->continuum.oConsumptionBg);
            conf.lookupValue("o-border", this->continuum.oBorder);

            conf.lookupValue("hCoarse", this->continuum.hCoarse);
            conf.lookupValue("hRefined", this->continuum.hRefined);
            conf.lookupValue("deltaT", this->continuum.deltaT);
        }
        catch(const SettingNotFoundException &nfex)
        {
            this->created_ = false;
            std::cerr << "Continuum Settings Not Found!" << std::endl;
        }
    }

    void getAgentValues()
    {
        const Setting& conf = this->cfg.getRoot()["agents"]["agent"];
        try
        {
            conf.lookupValue("nucleus-radius", this->agent.nucleusRadius);
            conf.lookupValue("radius", this->agent.radius);
            conf.lookupValue("action-radius", this->agent.actionRadius);
            conf.lookupValue("o-consumption", this->agent.oConsumption);
            this->agent.actionRadius *= this->agent.radius;
        }
        catch(const SettingNotFoundException &nfex)
        {
            this->created_ = false;
            std::cerr << "Agent Settings Not Found!" << std::endl;
        }
    }

    void getForcesValues()
    {
        const Setting& conf = this->cfg.getRoot()["agents"]["forces"];
        try
        {
            conf.lookupValue("c_cca", this->forces.c_cca);
            conf.lookupValue("c_ccr", this->forces.c_ccr);
            conf.lookupValue("K", this->forces.K);
            conf.lookupValue("c_ct", this->forces.c_ct);
            conf.lookupValue("c_rct", this->forces.c_rct);
            this->forces.c_ct *= this->forces.K;
            this->forces.c_rct *= this->forces.K;

        }
        catch(const SettingNotFoundException &nfex)
        {
            this->created_ = false;
            std::cerr << "Forces Settings Not Found!" << std::endl;
        }
    }

    void getParametersValues()
    {
        const Setting& conf = this->cfg.getRoot()["agents"]["parameters"];
        try
        {
            conf.lookupValue("tauA", this->parameters.tauA);
            conf.lookupValue("tauP", this->parameters.tauP);
            conf.lookupValue("tauG1", this->parameters.tauG1);
            conf.lookupValue("delta_tt", this->parameters.delta_tt);
            conf.lookupValue("alphaP", this->parameters.alphaP );
            conf.lookupValue("sigmaH", this->parameters.sigmaH);
            conf.lookupValue("alphaA", this->parameters.alphaA);
            conf.lookupValue("gammaA", this->parameters.gammaA);
        }
        catch(const SettingNotFoundException &nfex)
        {
            this->created_ = false;
            std::cerr << "Parameters Settings Not Found!" << std::endl;
        }
    }

public:
    ConfigHandler (std::string configFile = "config.cfg")
    {
        this->created_ = false;
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
        }

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
