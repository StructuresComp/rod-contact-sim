#ifndef SETINPUT_H
#define SETINPUT_H

#include <iostream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

#include "Option.h"
#include "eigenIncludes.h"

class setInput {
public:

    typedef std::map <std::string, Option> OptionMap;
    OptionMap m_options;

    setInput();

    ~setInput();

    template<typename T>
    int AddOption(const std::string &name, const std::string &desc, const T &def);

    Option *GetOption(const std::string &name);

    bool &GetBoolOpt(const std::string &name);

    int &GetIntOpt(const std::string &name);

    double &GetScalarOpt(const std::string &name);

    Vector3d &GetVecOpt(const std::string &name);

    string &GetStringOpt(const std::string &name);

    int LoadOptions(const char *filename);

    int LoadOptions(const std::string &filename) {
        return LoadOptions(filename.c_str());
    }

    int LoadOptions(int argc, char **argv);

private:
    double helixRadius;
    double helixPitch;
    double rodRadius;
    double youngM;
    double Poisson;
    double shearM;
    double deltaTime;
    double totalTime;
    double tol, stol;
    int ipc;

    double axisLengthInput;
    double deltaLengthInput;

    double distance;

    int numFlagella;
    int maxIter; // maximum number of iterations
    double density;
    double epsilon;
    Vector3d gVector;
    double viscosity;
    bool render;
    bool saveData;
    double omega;

    double mu;
    double nu;
    double delta;
    double col_limit;
    int line_search;
};

#include "setInput.tcc"

#endif // PROBLEMBASE_H
