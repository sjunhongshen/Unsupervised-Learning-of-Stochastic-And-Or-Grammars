#include <cstring>
#include <cstdio>
#include <cassert>
#include <iostream>  
#include <unordered_map>
#include <map>
#include <algorithm>

#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>
#include <fstream>

#include "utils.h"
#include "attributes.h"

using namespace cv;
using namespace std;
using namespace AOG_LIB;

//output the sample result to a file
vector<vector<pair<string, vector<double> > > > Sampling(string filename, const AOG<string, double>& aog,
														 VertexId sample_root, int n_samples, int n_gibbs)
{
	vector<vector<pair<string, vector<double> > > > samplingResult;
	std::ofstream ofs(filename);
    boost::archive::binary_oarchive oa(ofs);
	vector<pair<unordered_map<pair<Symbolic_State<string>, int>, vector<double>, AOG_LIB::pair_hash>,
				unordered_map<pair<Symbolic_State<string>, int>, vector<pair<Symbolic_State<string>, int> >, AOG_LIB::pair_hash> > > sampled_data;
	for(int i = 0; i < n_samples; i++)
	{
		vector<pair<VertexId, int> > res;
		double prob = 1;
		vector<unordered_map<pair<VertexId, int>, vector<double>, AOG_LIB::pair_hash> > configurations;
		if(i % (n_samples / 10) == 0)
			cout << "[" << i << "/" << n_samples << "] " << endl;
		shared_ptr<unordered_map<pair<VertexId, int>, vector<pair<VertexId, int>>, AOG_LIB::pair_hash>> pg =
			aog.Sample(sample_root, res, configurations, prob, n_gibbs);
		
		vector<pair<string, vector<double> > > faceSeq;
		for(auto re: res)
		{
			faceSeq.push_back(make_pair(aog.GetVertexContent(re.first)->GetState().GetContent(), configurations.back().at(re)));;
		}

		samplingResult.push_back(faceSeq);

		unordered_map<pair<Symbolic_State<string>, int>, vector<double>, AOG_LIB::pair_hash> state2attr;
		unordered_map<pair<Symbolic_State<string>, int>, vector<pair<Symbolic_State<string>, int> >, AOG_LIB::pair_hash> state2rules;
		//write to the file
		for(auto it = pg->cbegin(); it != pg->cend(); it++)
		{
			pair<Symbolic_State<string>, int> parent_pair = make_pair(aog.GetVertexContent(it->first.first)->GetState(), it->first.second);
			if (configurations.back().count(it->first))
				state2attr[parent_pair] = configurations.back().at(it->first);
			vector<pair<Symbolic_State<string>, int> > children;
			for(auto vid: it->second)
			{
				pair<Symbolic_State<string>, int> child_pair = make_pair(aog.GetVertexContent(vid.first)->GetState(), vid.second);
				if (configurations.back().count(vid))
					state2attr[child_pair] = configurations.back().at(vid);
				children.push_back(child_pair);
			}
			state2rules[parent_pair] = children;
		}
		sampled_data.push_back(make_pair(state2attr, state2rules));
	}
	oa << sampled_data;
	// FaceShapeType EyeType EyeType NoseType MouthType EarType EarType
	return samplingResult;
}


AOG<string, double> CreateFaceAOG()
{
    vector<Symbolic_Rule<string> > rules;

    Symbolic_State<string> root("Face", false);
	
	/////////////////////////////////////////////////////////////////
	// Face Composition
	/////////////////////////////////////////////////////////////////
    Symbolic_State<string> faceShape("FaceShape", false);
    Symbolic_State<string> eye("Eye", false);
    Symbolic_State<string> nose("Nose", false);
    Symbolic_State<string> mouth("Mouth", false);
    Symbolic_State<string> ear("Ear", false);
    vector<Symbolic_State<string> > faceComposition = {
		faceShape, eye, nose, mouth, ear
	};
    Symbolic_Rule<string> faceAndRule(root, faceComposition); 
    rules.push_back(faceAndRule);

	/////////////////////////////////////////////////////////////////
	// Eye Composition
	/////////////////////////////////////////////////////////////////
	Symbolic_State<string> leftEye("LeftEye", false);
	Symbolic_State<string> rightEye("RightEye", false);
	vector<Symbolic_State<string> > eyeComposition = {leftEye, rightEye};
	Symbolic_Rule<string> eyeAndRule(eye, eyeComposition); 
    rules.push_back(eyeAndRule);

	Symbolic_State<string> leftEyeType1("LEyeType1", true);
	Symbolic_State<string> leftEyeType2("LEyeType2", true);
	Symbolic_State<string> leftEyeType3("LEyeType3", true);
	Symbolic_State<string> leftEyeType4("LEyeType4", true);
	vector<Symbolic_State<string> > leftEyeBranch1 = {leftEyeType1};
	vector<Symbolic_State<string> > leftEyeBranch2 = {leftEyeType2};
	vector<Symbolic_State<string> > leftEyeBranch3 = {leftEyeType3};
	vector<Symbolic_State<string> > leftEyeBranch4 = {leftEyeType4};
	Symbolic_Rule<string> leftEyeOrRule1(leftEye, leftEyeBranch1);
	Symbolic_Rule<string> leftEyeOrRule2(leftEye, leftEyeBranch2);
	Symbolic_Rule<string> leftEyeOrRule3(leftEye, leftEyeBranch3);
	Symbolic_Rule<string> leftEyeOrRule4(leftEye, leftEyeBranch4);
	rules.push_back(leftEyeOrRule1);
	rules.push_back(leftEyeOrRule2);
	rules.push_back(leftEyeOrRule3);
	rules.push_back(leftEyeOrRule4);

	Symbolic_State<string> rightEyeType1("REyeType1", true);
	Symbolic_State<string> rightEyeType2("REyeType2", true);
	Symbolic_State<string> rightEyeType3("REyeType3", true);
	Symbolic_State<string> rightEyeType4("REyeType4", true);
	vector<Symbolic_State<string> > rightEyeBranch1 = {rightEyeType1};
	vector<Symbolic_State<string> > rightEyeBranch2 = {rightEyeType2};
	vector<Symbolic_State<string> > rightEyeBranch3 = {rightEyeType3};
	vector<Symbolic_State<string> > rightEyeBranch4 = {rightEyeType4};
	Symbolic_Rule<string> rightEyeOrRule1(rightEye, rightEyeBranch1);
	Symbolic_Rule<string> rightEyeOrRule2(rightEye, rightEyeBranch2);
	Symbolic_Rule<string> rightEyeOrRule3(rightEye, rightEyeBranch3);
	Symbolic_Rule<string> rightEyeOrRule4(rightEye, rightEyeBranch4);
	rules.push_back(rightEyeOrRule1);
	rules.push_back(rightEyeOrRule2);
	rules.push_back(rightEyeOrRule3);
	rules.push_back(rightEyeOrRule4);

	/////////////////////////////////////////////////////////////////
	// Ear Composition
	/////////////////////////////////////////////////////////////////
	Symbolic_State<string> leftEar("LeftEar", false);
	Symbolic_State<string> rightEar("RightEar", false);
	vector<Symbolic_State<string> > earComposition = {leftEar, rightEar};
	Symbolic_Rule<string> earAndRule(ear, earComposition); 
    rules.push_back(earAndRule);

	Symbolic_State<string> leftEarType1("LEarType1", true);
	Symbolic_State<string> leftEarType2("LEarType2", true);
	Symbolic_State<string> leftEarType3("LEarType3", true);
	vector<Symbolic_State<string> > leftEarBranch1 = {leftEarType1};
	vector<Symbolic_State<string> > leftEarBranch2 = {leftEarType2};
	vector<Symbolic_State<string> > leftEarBranch3 = {leftEarType3};
	Symbolic_Rule<string> leftEarOrRule1(leftEar, leftEarBranch1);
	Symbolic_Rule<string> leftEarOrRule2(leftEar, leftEarBranch2);
	Symbolic_Rule<string> leftEarOrRule3(leftEar, leftEarBranch3);
	rules.push_back(leftEarOrRule1);
	rules.push_back(leftEarOrRule2);
	rules.push_back(leftEarOrRule3);

	Symbolic_State<string> rightEarType1("REarType1", true);
	Symbolic_State<string> rightEarType2("REarType2", true);
	Symbolic_State<string> rightEarType3("REarType3", true);
	vector<Symbolic_State<string> > rightEarBranch1 = {rightEarType1};
	vector<Symbolic_State<string> > rightEarBranch2 = {rightEarType2};
	vector<Symbolic_State<string> > rightEarBranch3 = {rightEarType3};
	Symbolic_Rule<string> rightEarOrRule1(rightEar, rightEarBranch1);
	Symbolic_Rule<string> rightEarOrRule2(rightEar, rightEarBranch2);
	Symbolic_Rule<string> rightEarOrRule3(rightEar, rightEarBranch3);
	rules.push_back(rightEarOrRule1);
	rules.push_back(rightEarOrRule2);
	rules.push_back(rightEarOrRule3);
	
	/////////////////////////////////////////////////////////////////
	// Nose Composition
	/////////////////////////////////////////////////////////////////
	Symbolic_State<string> noseType1("NoseType1", true);
	Symbolic_State<string> noseType2("NoseType2", true);
	Symbolic_State<string> noseType3("NoseType3", true);
	vector<Symbolic_State<string> > noseBranch1 = {noseType1};
	vector<Symbolic_State<string> > noseBranch2 = {noseType2};
	vector<Symbolic_State<string> > noseBranch3 = {noseType3};
	Symbolic_Rule<string> noseOrRule1(nose, noseBranch1);
	Symbolic_Rule<string> noseOrRule2(nose, noseBranch2);
	Symbolic_Rule<string> noseOrRule3(nose, noseBranch3);
	rules.push_back(noseOrRule1);
	rules.push_back(noseOrRule2);
	rules.push_back(noseOrRule3);

	/////////////////////////////////////////////////////////////////
	// Mouth Composition
	/////////////////////////////////////////////////////////////////
	Symbolic_State<string> mouthType1("MouthType1", true);
	Symbolic_State<string> mouthType2("MouthType2", true);
	Symbolic_State<string> mouthType3("MouthType3", true);
	vector<Symbolic_State<string> > mouthBranch1 = {mouthType1};
	vector<Symbolic_State<string> > mouthBranch2 = {mouthType2};
	vector<Symbolic_State<string> > mouthBranch3 = {mouthType3};
	Symbolic_Rule<string> mouthOrRule1(mouth, mouthBranch1);
	Symbolic_Rule<string> mouthOrRule2(mouth, mouthBranch2);
	Symbolic_Rule<string> mouthOrRule3(mouth, mouthBranch3);
	rules.push_back(mouthOrRule1);
	rules.push_back(mouthOrRule2);
	rules.push_back(mouthOrRule3);

	/////////////////////////////////////////////////////////////////
	// Face Shape Composition
	/////////////////////////////////////////////////////////////////
	Symbolic_State<string> faceShapeType1("FaceShapeType1", true);
	Symbolic_State<string> faceShapeType2("FaceShapeType2", true);
	vector<Symbolic_State<string> > faceShapeBranch1 = {faceShapeType1};
	vector<Symbolic_State<string> > faceShapeBranch2 = {faceShapeType2};
	Symbolic_Rule<string> faceShapeOrRule1(faceShape, faceShapeBranch1);
	Symbolic_Rule<string> faceShapeOrRule2(faceShape, faceShapeBranch2);
	rules.push_back(faceShapeOrRule1);
	rules.push_back(faceShapeOrRule2);


	/////////////////////////////////////////////////////////////////
	// Create AOG
	/////////////////////////////////////////////////////////////////
    AOG<string, double> aog(rules);
    aog.SetRoot(root);
    
	/////////////////////////////////////////////////////////////////
	// Define Attributes and potential functions
	/////////////////////////////////////////////////////////////////
	GeometricProps* face_prop = new 
		GeometricProps (face_x_mean_idx, face_y_mean_idx, face_x_var, face_y_var, face_scale_var, face_pose_var);
	GeometricProps* nose_prop = new 
		GeometricProps (nose_x_mean_idx, nose_y_mean_idx, nose_x_var, nose_y_var, nose_scale_var, nose_pose_var);
	GeometricProps* mouth_prop = new 
		GeometricProps (mouth_x_mean_idx, mouth_y_mean_idx, mouth_x_var, mouth_y_var, mouth_scale_var, mouth_pose_var);
	GeometricProps* ear_left_prop = new 
		GeometricProps (ear_left_x_mean_idx, ear_left_y_mean_idx, ear_left_x_var, ear_left_y_var, ear_left_scale_var, ear_left_pose_var);
	GeometricProps* ear_right_prop = new 
		GeometricProps (ear_right_x_mean_idx, ear_right_y_mean_idx, ear_right_x_var, ear_right_y_var, ear_right_scale_var, ear_right_pose_var);
	GeometricProps* eye_left_prop = new 
		GeometricProps (eye_left_x_mean_idx, eye_left_y_mean_idx, eye_left_x_var, eye_left_y_var, eye_left_scale_var, eye_left_pose_var);
	GeometricProps* eye_right_prop = new 
		GeometricProps (eye_right_x_mean_idx, eye_right_y_mean_idx, eye_right_x_var, eye_right_y_var, eye_right_scale_var, eye_right_pose_var);

	vector<Symbolic_State<string> > leaf_states = aog.GetLeafStates();
	for(auto leaf_state: leaf_states)
    {
		aog.SetVertexAttributeFunc(aog.GetVertexIdByState(leaf_state), LeafAttributeFunc);
		aog.SetUnaryPotentialFunc(aog.GetVertexIdByState(leaf_state), LeafPotentialFunc);
		if (leaf_state.GetContent()[0] == 'F')
			SetupAttributesRange(aog, leaf_state, face_prop);
		else if (leaf_state.GetContent()[0] == 'N')
			SetupAttributesRange(aog, leaf_state, nose_prop);
		else if (leaf_state.GetContent()[0] == 'M')
			SetupAttributesRange(aog, leaf_state, mouth_prop);
		else if (leaf_state.GetContent().compare(0, 3, "LEa") == 0)
			SetupAttributesRange(aog, leaf_state, ear_left_prop);
		else if (leaf_state.GetContent().compare(0, 3, "REa") == 0)
			SetupAttributesRange(aog, leaf_state, ear_right_prop);
		else if (leaf_state.GetContent().compare(0, 3, "LEy") == 0)
			SetupAttributesRange(aog, leaf_state, eye_left_prop);
		else if (leaf_state.GetContent().compare(0, 3, "REy") == 0)
			SetupAttributesRange(aog, leaf_state, eye_right_prop);
	}
	aog.SetVertexAttributeFunc(aog.GetVertexIdByState(ear), SumAttributeFunc);
	aog.SetVertexAttributeFunc(aog.GetVertexIdByState(leftEar), SumAttributeFunc);
	aog.SetVertexAttributeFunc(aog.GetVertexIdByState(rightEar), SumAttributeFunc);
	aog.SetVertexAttributeFunc(aog.GetVertexIdByState(eye), SumAttributeFunc);
	aog.SetVertexAttributeFunc(aog.GetVertexIdByState(leftEye), SumAttributeFunc);
	aog.SetVertexAttributeFunc(aog.GetVertexIdByState(rightEye), SumAttributeFunc);
	aog.SetVertexAttributeFunc(aog.GetVertexIdByState(nose), SumAttributeFunc);
	aog.SetVertexAttributeFunc(aog.GetVertexIdByState(mouth), SumAttributeFunc);
	aog.SetVertexAttributeFunc(aog.GetVertexIdByState(faceShape), SumAttributeFunc);

	aog.SetBinaryPotentialFunc(make_pair(aog.GetVertexIdByState(leftEye), aog.GetVertexIdByState(faceShape)), RelativeToFacePotentialFunc);
	aog.SetBinaryPotentialParam(make_pair(aog.GetVertexIdByState(leftEye), aog.GetVertexIdByState(faceShape)), eye_left_prop);
	
	aog.SetBinaryPotentialFunc(make_pair(aog.GetVertexIdByState(rightEye), aog.GetVertexIdByState(faceShape)), RelativeToFacePotentialFunc);
	aog.SetBinaryPotentialParam(make_pair(aog.GetVertexIdByState(rightEye), aog.GetVertexIdByState(faceShape)), eye_right_prop);
	
	aog.SetBinaryPotentialFunc(make_pair(aog.GetVertexIdByState(leftEar), aog.GetVertexIdByState(faceShape)), RelativeToFacePotentialFunc);
	aog.SetBinaryPotentialParam(make_pair(aog.GetVertexIdByState(leftEar), aog.GetVertexIdByState(faceShape)), ear_left_prop);
	
	aog.SetBinaryPotentialFunc(make_pair(aog.GetVertexIdByState(rightEar), aog.GetVertexIdByState(faceShape)), RelativeToFacePotentialFunc);
	aog.SetBinaryPotentialParam(make_pair(aog.GetVertexIdByState(rightEar), aog.GetVertexIdByState(faceShape)), ear_right_prop);
	
	aog.SetBinaryPotentialFunc(make_pair(aog.GetVertexIdByState(nose), aog.GetVertexIdByState(faceShape)), RelativeToFacePotentialFunc);
	aog.SetBinaryPotentialParam(make_pair(aog.GetVertexIdByState(nose), aog.GetVertexIdByState(faceShape)), nose_prop);
	
	aog.SetBinaryPotentialFunc(make_pair(aog.GetVertexIdByState(mouth), aog.GetVertexIdByState(faceShape)), RelativeToFacePotentialFunc);
	aog.SetBinaryPotentialParam(make_pair(aog.GetVertexIdByState(mouth), aog.GetVertexIdByState(faceShape)), mouth_prop);
	/////////////////////////////////////////////////////////////////
	// Define Or-nodes weights
	/////////////////////////////////////////////////////////////////
    // Left Eye's weights
    unordered_map<VertexId, double> weights = aog.GetOutEdgeWeights(aog.GetVertexIdByState(leftEye), true);
    for (auto it = weights.begin(); it != weights.end(); it++)
    {
        if (aog.GetStateByVertexId(aog.ChildrenVertices(it->first)[0]) == leftEyeType1)
            it->second = 0.3;
        else if (aog.GetStateByVertexId(aog.ChildrenVertices(it->first)[0]) == leftEyeType2)
            it->second = 0.2;
        else if (aog.GetStateByVertexId(aog.ChildrenVertices(it->first)[0]) == leftEyeType3)
            it->second = 0.25;
		else
            it->second = 0.25;
    }
    aog.SetOutEdgeWeights(aog.GetVertexIdByState(leftEye), weights);

	// Right Eye's weights
	weights = aog.GetOutEdgeWeights(aog.GetVertexIdByState(rightEye), true);
    for (auto it = weights.begin(); it != weights.end(); it++)
    {
        if (aog.GetStateByVertexId(aog.ChildrenVertices(it->first)[0]) == rightEyeType1)
            it->second = 0.35;
        else if (aog.GetStateByVertexId(aog.ChildrenVertices(it->first)[0]) == rightEyeType2)
            it->second = 0.2;
        else if (aog.GetStateByVertexId(aog.ChildrenVertices(it->first)[0]) == rightEyeType3)
            it->second = 0.25;
		else
            it->second = 0.2;
    }
    aog.SetOutEdgeWeights(aog.GetVertexIdByState(rightEye), weights);

	// Left ear's weights
	weights = aog.GetOutEdgeWeights(aog.GetVertexIdByState(leftEar), true);
    for (auto it = weights.begin(); it != weights.end(); it++)
    {
        if (aog.GetStateByVertexId(aog.ChildrenVertices(it->first)[0]) == leftEarType1)
            it->second = 0.4;
        else if (aog.GetStateByVertexId(aog.ChildrenVertices(it->first)[0]) == leftEarType2)
            it->second = 0.2;
		else
            it->second = 0.4;
    }
    aog.SetOutEdgeWeights(aog.GetVertexIdByState(leftEar), weights);

	// Right ear's weights
	weights = aog.GetOutEdgeWeights(aog.GetVertexIdByState(rightEar), true);
    for (auto it = weights.begin(); it != weights.end(); it++)
    {
        if (aog.GetStateByVertexId(aog.ChildrenVertices(it->first)[0]) == rightEarType1)
            it->second = 0.4;
        else if (aog.GetStateByVertexId(aog.ChildrenVertices(it->first)[0]) == rightEarType2)
            it->second = 0.3;
		else
            it->second = 0.3;
    }
    aog.SetOutEdgeWeights(aog.GetVertexIdByState(rightEar), weights);

	// Nose's weights
	weights = aog.GetOutEdgeWeights(aog.GetVertexIdByState(nose), true);
    for (auto it = weights.begin(); it != weights.end(); it++)
    {
        if (aog.GetStateByVertexId(aog.ChildrenVertices(it->first)[0]) == noseType1)
            it->second = 0.35;
        else if (aog.GetStateByVertexId(aog.ChildrenVertices(it->first)[0]) == noseType2)
            it->second = 0.45;
		else
            it->second = 0.2;
    }
    aog.SetOutEdgeWeights(aog.GetVertexIdByState(nose), weights);

	// Mouth's weights
	weights = aog.GetOutEdgeWeights(aog.GetVertexIdByState(mouth), true);
    for (auto it = weights.begin(); it != weights.end(); it++)
    {
        if (aog.GetStateByVertexId(aog.ChildrenVertices(it->first)[0]) == mouthType1)
            it->second = 0.3;
        else if (aog.GetStateByVertexId(aog.ChildrenVertices(it->first)[0]) == mouthType2)
            it->second = 0.4;
		else
            it->second = 0.3;
    }
    aog.SetOutEdgeWeights(aog.GetVertexIdByState(mouth), weights);

	// FaceShape's weights
	weights = aog.GetOutEdgeWeights(aog.GetVertexIdByState(faceShape), true);
    for (auto it = weights.begin(); it != weights.end(); it++)
    {
        if (aog.GetStateByVertexId(aog.ChildrenVertices(it->first)[0]) == faceShapeType1)
            it->second = 0.7;
        else
            it->second = 0.3;
    }
    aog.SetOutEdgeWeights(aog.GetVertexIdByState(faceShape), weights);

	return aog;
}

void MkDir(const string &out)
{
	string cmd = "mkdir -p " + out;
	system(cmd.c_str());
}

int main(int argc, char *argv[])
{
	auto args = ArgParser(argc, argv);

	AOG<string, double> faceAOG = CreateFaceAOG();

	faceAOG.Visualize("./", "./cartoon_face_grammar_vis.txt");
	faceAOG.SaveGraph("./", "cartoon_face_grammar.txt");

	bool isNeg = stoi(args["--neg"]) ? true : false;
	string directory = args["-o"].length()? args["-o"]: "faces_gibbs:" + args["-g"];

	printf("Sampling %s %s samples at \\%s\n", args["-n"].c_str(), 
		isNeg ? "negative" : "positive", directory.c_str());
	
	MkDir(directory);
	vector<vector<pair<string, vector<double> > > > faces =
		Sampling(directory + "/face_details", faceAOG, faceAOG.GetRoot(), stoi(args["-n"]), stoi(args["-g"]));
	
	DrawFaces(directory, faces, isNeg);
	return 0;
}