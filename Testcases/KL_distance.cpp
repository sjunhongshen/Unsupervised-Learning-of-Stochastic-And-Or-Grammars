#include <string>
#include <iostream>
#include <unordered_map>
#include "AOG.h"

using namespace std;
using namespace AOG_LIB;

/*
	This function returns the distance of two AOGs.
	@params:
		aog1, aog2: two AOGs calcalated.
		calcType: the metric for distance
			type1(default): KL divergence : D(aog1 || aog2) = Sigma(P_aog1(x) * log(P_aog1(x) / P_aog2(x)))
			type2:
			type3:
		numOfSamples: num of Samples for each AOG
		epsilon: a small const to avoid zero division
	@returns:
		a double value dist
*/

template<class StateType, class AttributeType>
double CalculateDistance(AOG<StateType, AttributeType> aog1, AOG<StateType, AttributeType> aog2,
						 int calcType = 1, int numOfSamples = 100, double epsilon = 1e-9) {
	// Two maps used for counting
	unordered_map<StateType, int> map1;
	unordered_map<StateType, int> map2;
	vector<pair<VertexId, int> > res;
    vector<unordered_map<pair<VertexId, int>, vector<AttributeType>, AOG_LIB::pair_hash> > configurations;
	double prob = 1;
	//size of total counts
	int size1 = 0;
	int size2 = 0;
	for (int i = 0; i < numOfSamples; i++) {
		aog1.Sample(aog1.GetRoot(), res, configurations, prob);
		for (auto vid : res) {
			StateType state_content = aog1.GetStateByVertexId(vid.first).GetContent();
			if (map1.count(state_content) == 0) {
				map1[state_content] = 1;
			}
			else {
				map1[state_content] ++;
			}
			size1++;
		}

		aog2.Sample(aog2.GetRoot(), res, configurations, prob);
		for (auto vid : res) {
			StateType state_content = aog2.GetStateByVertexId(vid.first).GetContent();
			if (map2.count(state_content) == 0) {
				map2[state_content] = 1;
			}
			else {
				map2[state_content] ++;
			}
			size2++;
		}
	}
	if (calcType == 1) {
		double KL = 0;
		for (auto it = map1.begin(); it != map1.end(); it++) {
			StateType state_content = it->first;
			double prob1 = it->second / (double) size1;
			double prob2;
			if (map2.count(state_content) == 0) {
				prob2 = 0;
			}
			else {
				prob2 = map2[state_content] / (double) size2;
			}
			KL += prob1 * log(prob1 / (prob2 + epsilon));
		}
		return KL;
	}

	return 0;
}
AOG<string, string> createTestcase1() {
	vector<Symbolic_Rule<string> > rules;
	Symbolic_State<string> root("O", false);
	Symbolic_State<string> A("A", false);
	Symbolic_State<string> B("B", false);
	Symbolic_State<string> C("C", false);
	vector<Symbolic_State<string> > level1A = { A };
	vector<Symbolic_State<string> > level1B = { B };
	vector<Symbolic_State<string> > level1C = { C };
	Symbolic_Rule<string> level1Arule(root, level1A);
	Symbolic_Rule<string> level1Brule(root, level1B);
	Symbolic_Rule<string> level1Crule(root, level1C);
	rules.push_back(level1Arule);
	rules.push_back(level1Brule);
	rules.push_back(level1Crule);

	// A
	Symbolic_State<string> a1("a1", false);
	Symbolic_State<string> a2("a2", false);
	Symbolic_State<string> a3("a3", false);
	vector<Symbolic_State<string> > a1a2a3 = { a1, a2, a3 };
	Symbolic_Rule<string> A_a1a2a3(A, a1a2a3);
	rules.push_back(A_a1a2a3);

	Symbolic_State<string> a1A1("american", true);
	Symbolic_State<string> a1A2("engineer", true);
	Symbolic_State<string> a1B1("chinese", true);
	Symbolic_State<string> a1B2("student", true);
	vector<Symbolic_State<string> > a1A1a1A2 = { a1A1, a1A2 };
	vector<Symbolic_State<string> > a1B1a1B2 = { a1B1, a1B2 };
	Symbolic_Rule<string> a1_a1A1a1A2(a1, a1A1a1A2);
	Symbolic_Rule<string> a1_a1B1a1B2(a1, a1B1a1B2);
	rules.push_back(a1_a1A1a1A2);
	rules.push_back(a1_a1B1a1B2);

	Symbolic_State<string> a21("likes", true);
	Symbolic_State<string> a22("wants", true);
	vector<Symbolic_State<string> > va21 = { a21 };
	vector<Symbolic_State<string> > va22 = { a22 };
	Symbolic_Rule<string> a2_a21(a2, va21);
	Symbolic_Rule<string> a2_a22(a2, va22);

	rules.push_back(a2_a21);
	rules.push_back(a2_a22);


	Symbolic_State<string> a31("girl", true);
	Symbolic_State<string> a32("knowledge", true);
	vector<Symbolic_State<string> > va31 = { a31 };
	vector<Symbolic_State<string> > va32 = { a32 };
	Symbolic_Rule<string> a3_a31(a3, va31);
	Symbolic_Rule<string> a3_a32(a3, va32);
	rules.push_back(a3_a31);
	rules.push_back(a3_a32);


	// B
	Symbolic_State<string> b1("b1", false);
	Symbolic_State<string> b2("b2", false);
	vector<Symbolic_State<string> > b1b2 = { b1, b2 };
	Symbolic_Rule<string> B_b1b2(B, b1b2);
	rules.push_back(B_b1b2);

	Symbolic_State<string> B1A1("the", true);
	Symbolic_State<string> B1A2("green", true);
	Symbolic_State<string> B1A3("robot", true);
	Symbolic_State<string> B1B1("a", true);
	Symbolic_State<string> B1B2("cyborg", true);

	vector<Symbolic_State<string> > B1A1B1A2B1A3 = { B1A1, B1A2, B1A3 };
	vector<Symbolic_State<string> > B1B1B1B2 = { B1B1, B1B2 };
	Symbolic_Rule<string> b1_B1A1B1A2B1A3(b1, B1A1B1A2B1A3);
	Symbolic_Rule<string> b1_B1B1B1B2(b1, B1B1B1B2);
	rules.push_back(b1_B1A1B1A2B1A3);
	rules.push_back(b1_B1B1B1B2);

	Symbolic_State<string> B2A1("can", true);
	Symbolic_State<string> B2A2("dance", true);
	Symbolic_State<string> B2B1("never", true);
	Symbolic_State<string> B2B2("sleep", true);
	vector<Symbolic_State<string> > B2A1B2A2 = { B2A1, B2A2 };
	vector<Symbolic_State<string> > B2B1B2B2 = { B2B1, B2B2 };
	Symbolic_Rule<string> b2_B2A1B2A2(b2, B2A1B2A2);
	Symbolic_Rule<string> b2_B2B1B2B2(b2, B2B1B2B2);
	rules.push_back(b2_B2A1B2A2);
	rules.push_back(b2_B2B1B2B2);


	// C
	Symbolic_State<string> c1("c1", false);
	Symbolic_State<string> c2("c2", false);
	vector<Symbolic_State<string> > c1c2 = { c1, c2 };
	Symbolic_Rule<string> C_c1c2(C, c1c2);
	rules.push_back(C_c1c2);

	Symbolic_State<string> C1A1("Batman", true);
	Symbolic_State<string> C1A2("and", true);
	Symbolic_State<string> C1A3("Superman", true);
	Symbolic_State<string> C1B1("ninjia", true);
	Symbolic_State<string> C1B2("turtles", true);

	vector<Symbolic_State<string> > C1A1C1A2C1A3 = { C1A1, C1A2, C1A3 };
	vector<Symbolic_State<string> > C1B1C1B2 = { C1B1, C1B2 };
	Symbolic_Rule<string> c1_C1A1C1A2C1A3(c1, C1A1C1A2C1A3);
	Symbolic_Rule<string> c1_C1B1C1B2(c1, C1B1C1B2);
	rules.push_back(c1_C1A1C1A2C1A3);
	rules.push_back(c1_C1B1C1B2);

	Symbolic_State<string> C2A1("always", true);
	Symbolic_State<string> C2A2("fight", true);
	Symbolic_State<string> C2B1("live", true);
	Symbolic_State<string> C2B2("together", true);
	vector<Symbolic_State<string> > C2A1C2A2 = { C2A1, C2A2 };
	vector<Symbolic_State<string> > C2B1C2B2 = { C2B1, C2B2 };
	Symbolic_Rule<string> c2_C2A1C2A2(c2, C2A1C2A2);
	Symbolic_Rule<string> c2_C2B1C2B2(c2, C2B1C2B2);
	rules.push_back(c2_C2A1C2A2);
	rules.push_back(c2_C2B1C2B2);

	AOG<string, string> aog1(rules);
	aog1.SetRoot(root);
	return aog1;
}

AOG<string, string> createTestcase2() {
	vector<Symbolic_Rule<string> > rules;
	Symbolic_State<string> root("O", false);
	Symbolic_State<string> A("A", false);
	Symbolic_State<string> B("B", false);
	Symbolic_State<string> C("C", false);
	vector<Symbolic_State<string> > level1A = { A };
	vector<Symbolic_State<string> > level1B = { B };
	vector<Symbolic_State<string> > level1C = { C };
	Symbolic_Rule<string> level1Arule(root, level1A);
	Symbolic_Rule<string> level1Brule(root, level1B);
	Symbolic_Rule<string> level1Crule(root, level1C);
	rules.push_back(level1Arule);
	rules.push_back(level1Brule);
	rules.push_back(level1Crule);

	// A
	Symbolic_State<string> a1("a1", false);
	Symbolic_State<string> a2("a2", false);
	Symbolic_State<string> a3("a3", false);
	vector<Symbolic_State<string> > a1a2a3 = { a1, a2, a3 };
	Symbolic_Rule<string> A_a1a2a3(A, a1a2a3);
	rules.push_back(A_a1a2a3);

	Symbolic_State<string> a1A1("american", true);
	Symbolic_State<string> a1A2("chinese", true);
	Symbolic_State<string> a1B1("engineer", true);
	Symbolic_State<string> a1B2("student", true);
	vector<Symbolic_State<string> > a1A1a1A2 = { a1A1, a1A2 };
	vector<Symbolic_State<string> > a1B1a1B2 = { a1B1, a1B2 };
	Symbolic_Rule<string> a1_a1A1a1A2(a1, a1A1a1A2);
	Symbolic_Rule<string> a1_a1B1a1B2(a1, a1B1a1B2);
	rules.push_back(a1_a1A1a1A2);
	rules.push_back(a1_a1B1a1B2);

	Symbolic_State<string> a21("likes", true);
	Symbolic_State<string> a22("wants", true);
	vector<Symbolic_State<string> > va21 = { a21 };
	vector<Symbolic_State<string> > va22 = { a22 };
	Symbolic_Rule<string> a2_a21(a2, va21);
	Symbolic_Rule<string> a2_a22(a2, va22);

	rules.push_back(a2_a21);
	rules.push_back(a2_a22);


	Symbolic_State<string> a31("girl", true);
	Symbolic_State<string> a32("knowledge", true);
	vector<Symbolic_State<string> > va31 = { a31 };
	vector<Symbolic_State<string> > va32 = { a32 };
	Symbolic_Rule<string> a3_a31(a3, va31);
	Symbolic_Rule<string> a3_a32(a3, va32);
	rules.push_back(a3_a31);
	rules.push_back(a3_a32);


	// B
	Symbolic_State<string> b1("b1", false);
	Symbolic_State<string> b2("b2", false);
	vector<Symbolic_State<string> > b1b2 = { b1, b2 };
	Symbolic_Rule<string> B_b1b2(B, b1b2);
	rules.push_back(B_b1b2);

	Symbolic_State<string> B1A1("the", true);
	Symbolic_State<string> B1A2("robot", true);
	Symbolic_State<string> B1A3("green", true);
	Symbolic_State<string> B1B1("cyborg", true);
	Symbolic_State<string> B1B2("a", true);

	vector<Symbolic_State<string> > B1A1B1A2B1A3 = { B1A1, B1A2, B1A3 };
	vector<Symbolic_State<string> > B1B1B1B2 = { B1B1, B1B2 };
	Symbolic_Rule<string> b1_B1A1B1A2B1A3(b1, B1A1B1A2B1A3);
	Symbolic_Rule<string> b1_B1B1B1B2(b1, B1B1B1B2);
	rules.push_back(b1_B1A1B1A2B1A3);
	rules.push_back(b1_B1B1B1B2);

	Symbolic_State<string> B2A1("never", true);
	Symbolic_State<string> B2A2("dance", true);
	Symbolic_State<string> B2B1("can", true);
	Symbolic_State<string> B2B2("sleep", true);
	vector<Symbolic_State<string> > B2A1B2A2 = { B2A1, B2A2 };
	vector<Symbolic_State<string> > B2B1B2B2 = { B2B1, B2B2 };
	Symbolic_Rule<string> b2_B2A1B2A2(b2, B2A1B2A2);
	Symbolic_Rule<string> b2_B2B1B2B2(b2, B2B1B2B2);
	rules.push_back(b2_B2A1B2A2);
	rules.push_back(b2_B2B1B2B2);


	// C
	Symbolic_State<string> c1("c1", false);
	Symbolic_State<string> c2("c2", false);
	vector<Symbolic_State<string> > c1c2 = { c1, c2 };
	Symbolic_Rule<string> C_c1c2(C, c1c2);
	rules.push_back(C_c1c2);

	Symbolic_State<string> C1A1("Batman", true);
	Symbolic_State<string> C1A2("and", true);
	Symbolic_State<string> C1A3("Superman", true);
	Symbolic_State<string> C1B1("ninjia", true);
	Symbolic_State<string> C1B2("turtles", true);

	vector<Symbolic_State<string> > C1A1C1A2C1A3 = { C1A1, C1A2, C1A3 };
	vector<Symbolic_State<string> > C1B1C1B2 = { C1B1, C1B2 };
	Symbolic_Rule<string> c1_C1A1C1A2C1A3(c1, C1A1C1A2C1A3);
	Symbolic_Rule<string> c1_C1B1C1B2(c1, C1B1C1B2);
	rules.push_back(c1_C1A1C1A2C1A3);
	rules.push_back(c1_C1B1C1B2);

	Symbolic_State<string> C2A1("always", true);
	Symbolic_State<string> C2A2("fight", true);
	Symbolic_State<string> C2B1("live", true);
	Symbolic_State<string> C2B2("together", true);
	vector<Symbolic_State<string> > C2A1C2A2 = { C2A1, C2A2 };
	vector<Symbolic_State<string> > C2B1C2B2 = { C2B1, C2B2 };
	Symbolic_Rule<string> c2_C2A1C2A2(c2, C2A1C2A2);
	Symbolic_Rule<string> c2_C2B1C2B2(c2, C2B1C2B2);
	rules.push_back(c2_C2A1C2A2);
	rules.push_back(c2_C2B1C2B2);

	AOG<string, string> aog2(rules);
	aog2.SetRoot(root);
	return aog2;
}

AOG<string, string> createTestcase3() {

	vector<Symbolic_Rule<string> > rules;
	Symbolic_State<string> root("O", false);
	Symbolic_State<string> A("A", false);
	Symbolic_State<string> B("B", false);
	Symbolic_State<string> C("C", false);
	vector<Symbolic_State<string> > level1A = { A };
	vector<Symbolic_State<string> > level1B = { B };
	vector<Symbolic_State<string> > level1C = { C };
	Symbolic_Rule<string> level1Arule(root, level1A);
	Symbolic_Rule<string> level1Brule(root, level1B);
	Symbolic_Rule<string> level1Crule(root, level1C);
	rules.push_back(level1Arule);
	rules.push_back(level1Brule);
	rules.push_back(level1Crule);

	// A
	Symbolic_State<string> a1("a1", false);
	Symbolic_State<string> a2("a2", false);
	Symbolic_State<string> a3("a3", false);
	vector<Symbolic_State<string> > a1a2a3 = { a1, a2, a3 };
	Symbolic_Rule<string> A_a1a2a3(A, a1a2a3);
	rules.push_back(A_a1a2a3);

	Symbolic_State<string> a1A1("american!", true);
	Symbolic_State<string> a1A2("engineer!", true);
	Symbolic_State<string> a1B1("chinese!", true);
	Symbolic_State<string> a1B2("student!", true);
	vector<Symbolic_State<string> > a1A1a1A2 = { a1A1, a1A2 };
	vector<Symbolic_State<string> > a1B1a1B2 = { a1B1, a1B2 };
	Symbolic_Rule<string> a1_a1A1a1A2(a1, a1A1a1A2);
	Symbolic_Rule<string> a1_a1B1a1B2(a1, a1B1a1B2);
	rules.push_back(a1_a1A1a1A2);
	rules.push_back(a1_a1B1a1B2);

	Symbolic_State<string> a21("likes!", true);
	Symbolic_State<string> a22("wants!", true);
	vector<Symbolic_State<string> > va21 = { a21 };
	vector<Symbolic_State<string> > va22 = { a22 };
	Symbolic_Rule<string> a2_a21(a2, va21);
	Symbolic_Rule<string> a2_a22(a2, va22);

	rules.push_back(a2_a21);
	rules.push_back(a2_a22);


	Symbolic_State<string> a31("girl!", true);
	Symbolic_State<string> a32("knowledge!", true);
	vector<Symbolic_State<string> > va31 = { a31 };
	vector<Symbolic_State<string> > va32 = { a32 };
	Symbolic_Rule<string> a3_a31(a3, va31);
	Symbolic_Rule<string> a3_a32(a3, va32);
	rules.push_back(a3_a31);
	rules.push_back(a3_a32);


	// B
	Symbolic_State<string> b1("b1!", false);
	Symbolic_State<string> b2("b2!", false);
	vector<Symbolic_State<string> > b1b2 = { b1, b2 };
	Symbolic_Rule<string> B_b1b2(B, b1b2);
	rules.push_back(B_b1b2);

	Symbolic_State<string> B1A1("the!", true);
	Symbolic_State<string> B1A2("green!", true);
	Symbolic_State<string> B1A3("robot!", true);
	Symbolic_State<string> B1B1("a!", true);
	Symbolic_State<string> B1B2("cyborg!", true);

	vector<Symbolic_State<string> > B1A1B1A2B1A3 = { B1A1, B1A2, B1A3 };
	vector<Symbolic_State<string> > B1B1B1B2 = { B1B1, B1B2 };
	Symbolic_Rule<string> b1_B1A1B1A2B1A3(b1, B1A1B1A2B1A3);
	Symbolic_Rule<string> b1_B1B1B1B2(b1, B1B1B1B2);
	rules.push_back(b1_B1A1B1A2B1A3);
	rules.push_back(b1_B1B1B1B2);

	Symbolic_State<string> B2A1("can!", true);
	Symbolic_State<string> B2A2("dance!", true);
	Symbolic_State<string> B2B1("never!", true);
	Symbolic_State<string> B2B2("sleep!", true);
	vector<Symbolic_State<string> > B2A1B2A2 = { B2A1, B2A2 };
	vector<Symbolic_State<string> > B2B1B2B2 = { B2B1, B2B2 };
	Symbolic_Rule<string> b2_B2A1B2A2(b2, B2A1B2A2);
	Symbolic_Rule<string> b2_B2B1B2B2(b2, B2B1B2B2);
	rules.push_back(b2_B2A1B2A2);
	rules.push_back(b2_B2B1B2B2);


	// C
	Symbolic_State<string> c1("c1!", false);
	Symbolic_State<string> c2("c2!", false);
	vector<Symbolic_State<string> > c1c2 = { c1, c2 };
	Symbolic_Rule<string> C_c1c2(C, c1c2);
	rules.push_back(C_c1c2);

	Symbolic_State<string> C1A1("Batman!", true);
	Symbolic_State<string> C1A2("and!", true);
	Symbolic_State<string> C1A3("Superman!", true);
	Symbolic_State<string> C1B1("ninjia!", true);
	Symbolic_State<string> C1B2("turtles!", true);

	vector<Symbolic_State<string> > C1A1C1A2C1A3 = { C1A1, C1A2, C1A3 };
	vector<Symbolic_State<string> > C1B1C1B2 = { C1B1, C1B2 };
	Symbolic_Rule<string> c1_C1A1C1A2C1A3(c1, C1A1C1A2C1A3);
	Symbolic_Rule<string> c1_C1B1C1B2(c1, C1B1C1B2);
	rules.push_back(c1_C1A1C1A2C1A3);
	rules.push_back(c1_C1B1C1B2);

	Symbolic_State<string> C2A1("always!", true);
	Symbolic_State<string> C2A2("fight!", true);
	Symbolic_State<string> C2B1("live!", true);
	Symbolic_State<string> C2B2("together!", true);
	vector<Symbolic_State<string> > C2A1C2A2 = { C2A1, C2A2 };
	vector<Symbolic_State<string> > C2B1C2B2 = { C2B1, C2B2 };
	Symbolic_Rule<string> c2_C2A1C2A2(c2, C2A1C2A2);
	Symbolic_Rule<string> c2_C2B1C2B2(c2, C2B1C2B2);
	rules.push_back(c2_C2A1C2A2);
	rules.push_back(c2_C2B1C2B2);

	AOG<string, string> aog3(rules);
	aog3.SetRoot(root);
	return aog3;
}

int main()
{
	//Define an AOG
	AOG<string, string> aog1 = createTestcase1();
	AOG<string, string> aog2 = createTestcase2();
	AOG<string, string> aog3 = createTestcase3();
	// testcase 1: Same AOGs. Result: 2.06e-6
	cout << "The Distance for these two aogs are :" << CalculateDistance(aog1, aog1) << "\n";

	// testcase 2: Slightly change the state sequence. Result: 5.30e-3
	cout << "The Distance for these two aogs are :" << CalculateDistance(aog1, aog2) << "\n";

	// testcase 3: Completely different AOGs. Result: 17.4884
	cout << "The Distance for these two aogs are :" << CalculateDistance(aog1, aog3) << "\n";
}