#ifndef _UTILS_H
#define _UTILS_H

#include <cstdlib>
#include <map>

#include <opencv2/core/core.hpp>  
#include <opencv2/highgui/highgui.hpp>  
#include <opencv2/imgproc.hpp>

#include "components.h"
#include "argparse.hpp"

typedef std::vector<int> VI;
typedef std::vector<VI> VVI;

extern std::map<std::string, int> SymbolicType2Int;

void MkDir(const string &out);
void DrawFaces(const string &out, const vector<vector<pair<string, vector<double> > > >& faces, bool isNeg=false);
ap::argmap ArgParser(int argc, char *argv[]);

#endif 