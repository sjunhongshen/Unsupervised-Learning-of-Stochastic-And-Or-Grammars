#include "utils.h"

using namespace std;

map<string, int> SymbolicType2Int = {
	{"LEyeType1", 1},
	{"LEyeType2", 2},
	{"LEyeType3", 3},
	{"LEyeType4", 4},
	{"LEarType1", 1},
	{"LEarType2", 2},
	{"LEarType3", 3},
	{"REyeType1", 1},
	{"REyeType2", 2},
	{"REyeType3", 3},
	{"REyeType4", 4},
	{"REarType1", 1},
	{"REarType2", 2},
	{"REarType3", 3},
	{"NoseType1", 1},
	{"NoseType2", 2},
	{"NoseType3", 3},
	{"MouthType1", 1},
	{"MouthType2", 2},
	{"MouthType3", 3},
	{"FaceShapeType1", 1},
	{"FaceShapeType2", 2}
};

void DrawFaces(const string &out, const vector<vector<pair<string, vector<double> > > > &faces, bool isNeg)
{

	for(int i = 0; i < faces.size(); i++){
		// FaceShapeType EyeType EyeType NoseType MouthType EarType EarType
		Face face = Face(SymbolicType2Int.at(faces[i][0].first));
		face.copyXYSO(faces[i][0].second);
		Eye leftEye = Eye(SymbolicType2Int.at(faces[i][1].first), "Left");
		leftEye.copyXYSO(faces[i][1].second);
		Eye rightEye = Eye(SymbolicType2Int.at(faces[i][2].first), "Right");
		rightEye.copyXYSO(faces[i][2].second);
		Nose nose = Nose(SymbolicType2Int.at(faces[i][3].first));
		nose.copyXYSO(faces[i][3].second);
		Mouth mouth = Mouth(SymbolicType2Int.at(faces[i][4].first));
		mouth.copyXYSO(faces[i][4].second);
		Ear leftEar = Ear(SymbolicType2Int.at(faces[i][5].first), "Left");
		leftEar.copyXYSO(faces[i][5].second);
		Ear rightEar = Ear(SymbolicType2Int.at(faces[i][6].first), "Right");
		rightEar.copyXYSO(faces[i][6].second);

		if(isNeg){
			leftEye.setNegSampleFlag(true);
			rightEye.setNegSampleFlag(true);
			nose.setNegSampleFlag(true);
			mouth.setNegSampleFlag(true);
			leftEar.setNegSampleFlag(true);
			rightEar.setNegSampleFlag(true);
		}

		// draw on panel
		Mat panel = Mat(WINDOW_SIZE, WINDOW_SIZE, CV_8UC1, Scalar(255));
		panel = face.drawPic(panel);
		panel = leftEye.drawPic(panel);
		panel = rightEye.drawPic(panel);
		panel = nose.drawPic(panel);
		panel = mouth.drawPic(panel);
		panel = leftEar.drawPic(panel);
		panel = rightEar.drawPic(panel);

		std::string filename = out + "/face_" + to_string(i) + ".png";
		imwrite(filename, panel);
	}
}

ap::argmap ArgParser(int argc, char *argv[])
{
	ap::parser p(argc, argv);
	
	p.add("-b", "--neg", "Sampling negative samples", ap::mode::BOOLEAN);
	p.add("-o", "--out", "Output directory", ap::mode::OPTIONAL);
	p.add("-n", "--n_samples", "Number of samples", ap::mode::REQUIRED);
	p.add("-g", "--n_gibbs", "Number of samples", ap::mode::REQUIRED);

	auto args = p.parse();

	if (!args.parsed_successfully()) {
		std::cerr << "Unsuccessful parse\n";	
		exit(1);
    }

	return args;
}