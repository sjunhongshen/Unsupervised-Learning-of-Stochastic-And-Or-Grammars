//
// Created by Luyao Yuan on 18/4/8.
// Compile this file with options -lboost_system -lboost_filesystem

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <chrono>

#include "AOG.h"
#include "Earley_Parser.h"
#include "structure_learner.h"

using namespace std;
using namespace AOG_LIB;
template <class T>
using SequenceType = std::vector<AOG_LIB::Symbolic_State<T> >;

int SAMPLING = 500;

void OutputToFile(string filename, const AOG<string, double>& t_aog, VertexId sample_root)
{
    vector<Symbolic_Rule<string> > all_rules = t_aog.GetRules();
    Symbolic_State<string> root_state = t_aog.GetStateByVertexId(t_aog.GetRoot());
    //open file
    ofstream outFile;
    outFile.open(filename);   
    if(outFile.is_open())
    {
        while (SAMPLING)
        {
            vector<pair<VertexId, int> > res;
            double sample_prob = 1;
            vector<unordered_map<pair<VertexId, int>, vector<double>, AOG_LIB::pair_hash> > configurations;
            shared_ptr<unordered_map<pair<VertexId, int>, vector<pair<VertexId, int>>, AOG_LIB::pair_hash>> pg
                = t_aog.Sample(sample_root, res, configurations, sample_prob);
            SAMPLING--;
            vector<Symbolic_State<string> > seq;
            for(auto vid : res)
            {
                outFile << t_aog.GetStateByVertexId(vid.first).GetContent() << " ";
                seq.push_back(t_aog.GetStateByVertexId(vid.first));
            }
            outFile<<"\n";
        }
    }
    else
    {
        std::cerr<<"Error opening file!\n";
        throw exception();
    }
    outFile.close();
}


bool copyFile(const char *SRC, const char *DEST)
{
	std::ifstream src(SRC, std::ios::binary);
	std::ofstream dest(DEST, std::ios::binary);
	dest << src.rdbuf();
	return src && dest;
}

int main(int argc, char *argv[])
{
	
	if (argc < 2)
	{
		std::cout << "[Argument needed!]" << std::endl;
		assert(0);
	}

	string FILEPATH = argv[1];

	string filename(FILEPATH);
	vector<SequenceType<string>> all_dataset = FileParser<string>(filename);
    int train_size = all_dataset.size() * 0.7;
    vector<SequenceType<string>> dataset(all_dataset.begin(), all_dataset.begin() + train_size);
    vector<SequenceType<string>> testing_dataset(all_dataset.begin()+train_size, all_dataset.end());

    unordered_map<SequenceType<string>, unsigned> rules_count;
    Symbolic_State<string> root("O", false);
    vector<Symbolic_Rule<string>> rules;
    for (unsigned i = 0; i < dataset.size(); i++)
    {
        rules_count[dataset[i]]++;
        // cerr<<dataset[i][0].GetContent()<<endl;
        Symbolic_Rule<string> new_rule(root, dataset[i]);
        rules.push_back(new_rule);
    }
    int test_iter = 0;
    while(test_iter < testing_dataset.size())
    {
        if(rules_count.count(testing_dataset[test_iter]) > 0)
            testing_dataset.erase(testing_dataset.begin() + test_iter);
        else
            test_iter++;
    }
    cout << "Read in dataset with " << dataset.size() << " training sequences"
         <<" and " << testing_dataset.size() << " testing data\n";

    AOG<string, double> aog(rules);
    aog.SetRoot(root);

    //set up initial weights of the data
    unordered_map<VertexId ,double> weights = aog.GetOutEdgeWeights(aog.GetRoot(), false);
    for(auto it = weights.begin(); it != weights.end(); it++)
    {
        vector<VertexId> children_id = aog.ChildrenVertices(it->first);
        SequenceType<string> children;
        for(auto cid: children_id)
            children.push_back(aog.GetVertexContent(cid)->GetState());
        weights[it->first] = rules_count.at(children);
        // cerr << "weights "<<it->first<<" "<<weights[it->first]<<endl;
    }
    aog.SetOutEdgeWeights(aog.GetRoot(), weights);
    double alpha = 2;
    Learn(aog, rules_count, alpha, 50, 20);
	aog.TruncateGraph();

    //check recall
    vector<Symbolic_Rule<string> > all_rules = aog.GetRules();
    shared_ptr<AOG_LIB_UTIL::grammar<string> > g =
        make_shared<AOG_LIB_UTIL::grammar<string> >(all_rules, vector<Symbolic_State<string> >{root});
    shared_ptr<AOG_LIB_UTIL::EarleyParser<string, double> > parser = make_shared<AOG_LIB_UTIL::EarleyParser<string, double> >(*g);
    int total_success = 0;
    for(int i = 0; i < testing_dataset.size(); i++)
    {
        if(i % (testing_dataset.size() / 10) == 0)
            cout << "[Parsing] testing data: " << i << endl;
        auto test_data = testing_dataset[i];
        for (auto t:test_data)
            cerr <<t.GetContent()<<" ";
        cerr<<endl;
        bool parsing_success = false;
        int pos = parser->parse(test_data.begin(), test_data.end(), std::cout, parsing_success);
        cerr<<pos<<endl;
        total_success += (pos == test_data.size());
    }
    if(!testing_dataset.empty())
        printf("Testing Data Parsing Accuracy: %f\n", 1.0 * total_success / testing_dataset.size());
    OutputToFile("./learned_grammar_samples.txt", aog, aog.GetVertexIdByState(root));
    aog.Visualize("./", "./test_learner_vis.txt");
	return 0;
}