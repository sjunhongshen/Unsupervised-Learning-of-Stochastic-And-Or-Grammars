//
// Created by Hao Wang on 18/02/22.
//

#include <string>
#include <iostream>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>

#include "AOG.h"

using namespace std;
using namespace AOG_LIB;

template <class T>
using SequenceType = std::vector<AOG_LIB::Symbolic_State<T> >;

using SSType = Symbolic_State<string>;
using SRType = Symbolic_Rule<string>;

std::default_random_engine generator;

vector<double> SimpleAttributeFunc(const AOG_Vertex<string, double> &vertex_self,
                                   const vector<vector<double>> *neighbor_attrs,
                                   void *args);

vector<double> LeafAttributeFunc(const AOG_Vertex<string, double>& vertex_self,
                                 const vector<vector<double> >* neighbor_attrs,
                                 void* args)
{
    vector<vector<double>> ranges = vertex_self.GetAttributesRange();
    vector<double> new_attr;
    int* int_ptr = static_cast<int*>(args);
    int attribute_pos = int_ptr[0];
    // randomly initialize attributes of this leaf node
    if(attribute_pos == -1)
    {
        // assert(vertex_self.GetAttributesRange().size());
        for (auto range: ranges)
        {
            uniform_int_distribution<int> dist(0, range.size() - 1);
            int attr_idx = dist(generator);
            new_attr.push_back(range[attr_idx]);
        }
    }
    else
    {
        new_attr = (*neighbor_attrs)[0];
        new_attr[attribute_pos]= ranges[attribute_pos][int_ptr[1]];
    }
    
    return new_attr;
}

vector<double> RootAttributeFunc(const AOG_Vertex<string, double>& vertex_self,
                                 const vector<vector<double> >* neighbor_attrs,
                                 void* args)
{
    if (!neighbor_attrs->empty())
        return SimpleAttributeFunc(vertex_self, neighbor_attrs, args);
    else
    {
        vector<vector<double>> ranges = vertex_self.GetAttributesRange();
        vector<double> new_attr;
        // assert(vertex_self.GetAttributesRange().size());
        for (auto range: ranges)
        {
            uniform_int_distribution<int> dist(0, range.size() - 1);
            int attr_idx = dist(generator);
            new_attr.push_back(range[attr_idx]);
        }
        return new_attr;
    }
}

vector<double> SimpleAttributeFunc(const AOG_Vertex<string, double> &vertex_self,
                                   const vector<vector<double>> *neighbor_attrs,
                                   void *args)
{
    vector<double> new_attr;
    int num_attributes = (*neighbor_attrs)[0].size();
    for (int i = 0; i < num_attributes; i++)
    {
        double sum = 0;
        for (auto n_attr: *neighbor_attrs)
            sum += n_attr[i];
        sum /= (*neighbor_attrs).size();
        new_attr.push_back(sum);
    }
    return new_attr;
}

double LeafPotentialFunc(const vector<double>& self_attr, void* args)
{
    return 1000 * (self_attr[0] - self_attr[1]) * (self_attr[0] - self_attr[1]);
}

double SimplePotentialFunc(const vector<vector<double> >& attributes, void* args)
{
    return 0 * (attributes[0][0] + attributes[0][1] - attributes[1][0] - attributes[1][1]) *
               (attributes[0][0] + attributes[0][1] - attributes[1][0] - attributes[1][1]);
}

//output the sample result to a file
void SampleToFile(string filename, const AOG<string, double>& t_aog, VertexId sample_root, int num_samples = 100)
{
    //open file
    ofstream outFile;
    outFile.open(filename + ".txt");
    std::ofstream ofs(filename + ".dat");
    boost::archive::binary_oarchive oa(ofs);
	vector<pair<unordered_map<pair<Symbolic_State<string>, int>, vector<double>, AOG_LIB::pair_hash>,
				unordered_map<pair<Symbolic_State<string>, int>, vector<pair<Symbolic_State<string>, int> >, AOG_LIB::pair_hash> > > sampled_data;
    if(outFile.is_open())
    {
        while (num_samples)
        {
            vector<pair<VertexId, int> > res;
            double sample_prob = 1;
            vector<unordered_map<pair<VertexId, int>, vector<double>, AOG_LIB::pair_hash> > configurations;
            cout << "[" << num_samples << "]";
            shared_ptr<unordered_map<pair<VertexId, int>, vector<pair<VertexId, int>>, AOG_LIB::pair_hash>> pg
                = t_aog.Sample(sample_root, res, configurations, sample_prob);
            num_samples--;
            unsigned num_depth = pg->size();
            for(auto it = pg->cbegin(); it != pg->cend(); it++)
            {
                if (configurations.back().count(it->first))
                    outFile << "<(" << it->first.first << "-" << it->first.second << "):" << t_aog.GetStateByVertexId(it->first.first).GetContent() << ":"
                            << configurations.back().at(it->first) << "> -> ";
                else
                    outFile << "<(" << it->first.first << "-" << it->first.second << "):" << t_aog.GetStateByVertexId(it->first.first).GetContent() << "> -> ";
                for(auto vid: it->second)
                {
                    if (configurations.back().count(vid))
                        outFile << "<(" << vid.first << "-" << vid.second << "):" << t_aog.GetStateByVertexId(vid.first).GetContent()
                                << ":" << configurations.back().at(vid) << "> ";
                    else
                        outFile << "<(" << vid.first << "-" << vid.second << "):" << t_aog.GetStateByVertexId(vid.first).GetContent()
                                << ":> ";
                }
                outFile << '\n';
            }
            for (auto vid : res)
                outFile << t_aog.GetStateByVertexId(vid.first).GetContent() << " ";
            outFile << "\n";
            unordered_map<pair<Symbolic_State<string>, int>, vector<double>, AOG_LIB::pair_hash> state2attr;
            unordered_map<pair<Symbolic_State<string>, int>, vector<pair<Symbolic_State<string>, int> >, AOG_LIB::pair_hash> state2rules;
            //write to the file
            for(auto it = pg->cbegin(); it != pg->cend(); it++)
            {
                pair<Symbolic_State<string>, int> parent_pair = make_pair(t_aog.GetVertexContent(it->first.first)->GetState(), it->first.second);
                if (configurations.back().count(it->first))
                    state2attr[parent_pair] = configurations.back().at(it->first);
                vector<pair<Symbolic_State<string>, int> > children;
                for(auto vid: it->second)
                {
                    pair<Symbolic_State<string>, int> child_pair = make_pair(t_aog.GetVertexContent(vid.first)->GetState(), vid.second);
                    if (configurations.back().count(vid))
                        state2attr[child_pair] = configurations.back().at(vid);
                    children.push_back(child_pair);
                }
                state2rules[parent_pair] = children;
            }
            sampled_data.push_back(make_pair(state2attr, state2rules));
        }
        oa << sampled_data;
    }
    else
    {
        std::cerr<<"Error opening file!\n";
        throw exception();
    }
    ofs.close();
    outFile.close();
}

int main()
{
    //Define an AOG
    vector<SRType> rules;
    SSType root("O", false);
    SSType A("A", false);
    SSType B("B", false);
    SSType C("C", false);
    vector<SSType> level1A = {A};
    vector<SSType> level1B = {B};
    vector<SSType> level1C = {C};
    SRType level1Arule(root, level1A); 
    SRType level1Brule(root, level1B); 
    SRType level1Crule(root, level1C);
    rules.push_back(level1Arule);
    rules.push_back(level1Brule); 
    rules.push_back(level1Crule); 
    
    // A
    SSType a1("a1", false);
    SSType a1_("a1_", false);
    SSType a2("a2", false);
    SSType a3("a3", false);
    vector<SSType> a1a2a3 = {a1, a2, a3};
    SRType A_a1a2a3(A, a1a2a3);
    rules.push_back(A_a1a2a3);

    SSType a1A1("a1A1", true);
    SSType a1A2("a1A2", true);
    SSType a1B1("a1B1", true);
    SSType a1B2("a1B2", true);
    vector<SSType> a1A1a1A2 = {a1A1, a1A2};
    vector<SSType> a1B1a1B2 = {a1B1, a1B2};
    SRType a1__a1B1a1B2(a1_, a1B1a1B2);
    SRType a1_a1_(a1, {a1_});
    SRType a1_a1A1a1A2(a1, a1A1a1A2);
    rules.push_back(a1_a1A1a1A2);
    rules.push_back(a1__a1B1a1B2);
    rules.push_back(a1_a1_);
    
    SSType a21("a21", true);
    SSType a22("a22", true);
    vector<SSType> va21 = {a21};
    vector<SSType> va22 = {a22};
    SRType a2_a21(a2, va21);
    SRType a2_a22(a2, va22);
    SRType a2_a1_(a2, {a1_});
    
    rules.push_back(a2_a21);
    rules.push_back(a2_a22);
    rules.push_back(a2_a1_);
    

    SSType a31("a31", true);
    SSType a32("a32", true);
    vector<SSType> va31 = {a31};
    vector<SSType> va32 = {a32};
    SRType a3_a31(a3, va31);
    SRType a3_a32(a3, va32);
    rules.push_back(a3_a31);
    rules.push_back(a3_a32);
    

    // B
    SSType b1("b1", false);
    SSType b2("b2", false);
    vector<SSType> b1b2 = {b1, b2};
    SRType B_b1b2(B, b1b2);
    rules.push_back(B_b1b2);

    SSType B1A1("b1A1", true);
    SSType B1A2("b1A2", true);
    SSType B1A3("b1A3", true);
    SSType B1B1("b1B1", true);
    SSType B1B2("b1B2", true);
    
    vector<SSType> B1A1B1A2B1A3 = {B1A1, B1A2, B1A3};
    vector<SSType> B1B1B1B2 = {B1B1, B1B2};
    SRType b1_B1A1B1A2B1A3(b1, B1A1B1A2B1A3);
    SRType b1_B1B1B1B2(b1, B1B1B1B2);
    rules.push_back(b1_B1A1B1A2B1A3);
    rules.push_back(b1_B1B1B1B2);

    SSType B2A1("b2A1", true);
    SSType B2A2("b2A2", true);
    SSType B2B1("b2B1", true);
    SSType B2B2("b2B2", true);
    vector<SSType> B2A1B2A2 = {B2A1, B2A2};
    vector<SSType> B2B1B2B2 = {B2B1, B2B2};
    SRType b2_B2A1B2A2(b2, B2A1B2A2);
    SRType b2_B2B1B2B2(b2, B2B1B2B2);
    rules.push_back(b2_B2A1B2A2);
    rules.push_back(b2_B2B1B2B2);


    // C
    SSType c1("c1", false);
    SSType c2("c2", false);
    vector<SSType> c1c2 = {c1, c2};
    SRType C_c1c2(C, c1c2);
    rules.push_back(C_c1c2);

    SSType C1A1("c1A1", true);
    SSType C1A2("c1A2", true);
    SSType C1A3("c1A3", true);
    SSType C1B1("c1B1", true);
    SSType C1B2("c1B2", true);
    
    vector<SSType> C1A1C1A2C1A3 = {C1A1, C1A2, C1A3};
    vector<SSType> C1B1C1B2 = {C1B1, C1B2};
    SRType c1_C1A1C1A2C1A3(c1, C1A1C1A2C1A3);
    SRType c1_C1B1C1B2(c1, C1B1C1B2);
    rules.push_back(c1_C1A1C1A2C1A3);
    rules.push_back(c1_C1B1C1B2);

    SSType C2A1("c2A1", true);
    SSType C2A2("c2A2", true);
    SSType C2B1("c2B1", true);
    SSType C2B2("c2B2", true);
    vector<SSType> C2A1C2A2 = {C2A1, C2A2};
    vector<SSType> C2B1C2B2 = {C2B1, C2B2};
    SRType c2_C2A1C2A2(c2, C2A1C2A2);
    SRType c2_C2B1C2B2(c2, C2B1C2B2);
    rules.push_back(c2_C2A1C2A2);
    rules.push_back(c2_C2B1C2B2);

    AOG<string, double> t_aog(rules);
    t_aog.SetRoot(root);
    VertexId root_id = t_aog.GetVertexIdByState(root);

    cout << t_aog.GetVertexContent(t_aog.GetVertexIdByState(C2A1))->IsNotOr() << endl;

    vector<SSType> leaf_states = t_aog.GetLeafStates();
    vector<SSType> all_states = t_aog.GetStates();
    vector<VertexId> non_terminal_vertices;
    vector<VertexId> leaf_vertices;
    for(auto as: all_states)
        if (!as.GetIsBasic())
            non_terminal_vertices.push_back(t_aog.GetVertexIdByState(as));
    for(auto leaf_state: leaf_states)
        leaf_vertices.push_back(t_aog.GetVertexIdByState(leaf_state));

    vector<double> digits;
    double minD = 0;
    double maxD = 1;
    double step = 0.01;
    while(minD <= maxD)
    {
        digits.push_back(minD);
        minD += step;
    }

    for (auto vid: leaf_vertices)
    {
        t_aog.SetVertexAttributesRanges(vid, {digits, digits});
        t_aog.SetVertexAttributeFunc(vid, LeafAttributeFunc);
        t_aog.SetUnaryPotentialFunc(vid, LeafPotentialFunc);
    }
    
    for (auto vid: non_terminal_vertices)
        if (t_aog.GetStateByVertexId(vid).GetContent()[0] != 'c' &&
            t_aog.GetStateByVertexId(vid).GetContent()[0] != 'C')
            t_aog.SetVertexAttributeFunc(vid, SimpleAttributeFunc);
    
    t_aog.SetVertexAttributeFunc(root_id, RootAttributeFunc);
    t_aog.SetVertexAttributesRanges(root_id, {{digits[0]}, {digits[1], digits[2]}});
    t_aog.SetBinaryPotentialFunc(make_pair(t_aog.GetVertexIdByState(b1), t_aog.GetVertexIdByState(b2)), SimplePotentialFunc);

    // A's weight
    unordered_map<VertexId, double> weights = t_aog.GetOutEdgeWeights(t_aog.GetVertexIdByState(root), true);
    for (auto it = weights.begin(); it != weights.end(); it++)
    {
        if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == A)
            it->second = 0.3;
        if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == B)
            it->second = 0.5;
        if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == C)
            it->second = 0.2;
    }
    t_aog.SetOutEdgeWeights(t_aog.GetVertexIdByState(root), weights);


    weights = t_aog.GetOutEdgeWeights(t_aog.GetVertexIdByState(a1), true);
    for (auto it = weights.begin(); it != weights.end(); it++)
    {
        if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == a1A1)
            it->second = 0.7;
        if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == a1B1)
            it->second = 0.3;
    }
    t_aog.SetOutEdgeWeights(t_aog.GetVertexIdByState(a1), weights);   


    weights = t_aog.GetOutEdgeWeights(t_aog.GetVertexIdByState(a2), true);
    for (auto it = weights.begin(); it != weights.end(); it++)
    {
        if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == a21)
            it->second = 0.46;
        if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == a22)
            it->second = 0.54;
    }
    t_aog.SetOutEdgeWeights(t_aog.GetVertexIdByState(a2), weights);


    weights = t_aog.GetOutEdgeWeights(t_aog.GetVertexIdByState(a3), true);
    for (auto it = weights.begin(); it != weights.end(); it++)
    {
        if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == a31)
            it->second = 0.76;
        if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == a32)
            it->second = 0.24;
    }
    t_aog.SetOutEdgeWeights(t_aog.GetVertexIdByState(a3), weights);
    
    // B's weight
    weights = t_aog.GetOutEdgeWeights(t_aog.GetVertexIdByState(b1), true);
    for (auto it = weights.begin(); it != weights.end(); it++)
    {
        if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == B1A1)
            it->second = 0.21;
        if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == B1B1)
            it->second = 0.79;
    }
    t_aog.SetOutEdgeWeights(t_aog.GetVertexIdByState(b1), weights);

    weights = t_aog.GetOutEdgeWeights(t_aog.GetVertexIdByState(b2), true);
    for (auto it = weights.begin(); it != weights.end(); it++)
    {
        if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == B2A1)
            it->second = 0.64;
        if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == B2B1)
            it->second = 0.36;
    }
    t_aog.SetOutEdgeWeights(t_aog.GetVertexIdByState(b2), weights);


    // C's weight
    weights = t_aog.GetOutEdgeWeights(t_aog.GetVertexIdByState(c1), true);
    for (auto it = weights.begin(); it != weights.end(); it++)
    {
        if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == C1A1)
            it->second = 0.71;
        if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == C1B1)
            it->second = 0.29;
    }
    t_aog.SetOutEdgeWeights(t_aog.GetVertexIdByState(c1), weights);

    weights = t_aog.GetOutEdgeWeights(t_aog.GetVertexIdByState(c2), true);
    for (auto it = weights.begin(); it != weights.end(); it++)
    {
        if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == C2A1)
            it->second = 0.34;
        if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == C2B1)
            it->second = 0.66;
    }
    t_aog.SetOutEdgeWeights(t_aog.GetVertexIdByState(c2), weights);

    t_aog.SaveGraph("./", "grammar_example_exp1.txt");
    t_aog.Visualize("./", "./grammar_example_exp1_vis.txt");

    SampleToFile("./grammar_example_exp1_samples", t_aog, t_aog.GetVertexIdByState(root), 3);
    return 0;
}