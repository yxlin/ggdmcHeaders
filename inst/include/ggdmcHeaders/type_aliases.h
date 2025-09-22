#pragma once
#include <map>
#include <string>
#include <vector>

using bool3D = std::vector<std::vector<std::vector<bool>>>;
using bool2D = std::vector<std::vector<bool>>;
using bool1D = std::vector<bool>;

using uint2D = std::vector<std::vector<unsigned int>>;
using uint1D = std::vector<unsigned int>;

using strVec = std::vector<std::string>;
using strMap = std::map<std::string, std::string>;

using MapStrVec = std::map<std::string, strVec>;
using MapStrStr = std::map<std::string, std::map<std::string, std::string>>;
using MapStrDbl = std::map<std::string, double>;

using double2D = std::vector<std::vector<double>>;
