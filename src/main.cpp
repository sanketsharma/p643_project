#include <simulate.h>
#include <fstream>
#include <sstream>
#include <easylogging++.h>

INITIALIZE_EASYLOGGINGPP

std::map<std::string, std::string> getConfigValues(const std::string& fname)
{
    std::map<std::string, std::string> configMap;
    std::ifstream isFile(fname.c_str());
    std::string line;
    while(std::getline(isFile, line))
    {
        std::istringstream isLine(line);
        std::string key;
        if(std::getline(isLine >> std::ws, key, '='))
        {
            std::string value;
            if(std::getline(isLine >> std::ws, value))
            {
                configMap[key]=value;
            }
        }
    }
    return configMap;
}

int main()
{
    el::Configurations conf("../config/logging.cfg");
    auto configValues = getConfigValues("../config/params.cfg");
    LOG(INFO) << "Starting simulation" ;
    p643::simulate(configValues);
    LOG(INFO) << "Ending simulation" ;
    return 0;
}
