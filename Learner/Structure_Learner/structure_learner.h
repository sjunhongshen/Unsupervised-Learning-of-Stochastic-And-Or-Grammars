#pragma once

#include <stdexcept>
#include <boost/numeric/ublas/matrix.hpp>
#include <queue>
#include <vector>
#include <iterator>
#include <cmath>
#include <stdio.h>

#include "AOG.h"
#include "And_Or_Fragment.h"
#include "Context_Matrix.h"

using namespace std;
using namespace AOG_LIB;

//update current data with a selected aof and its reduction positions
template <class StateType>
void UpdateData(const AOFStruct<StateType>& aof,
                const unordered_map<SequenceType<StateType>, vector<unsigned>>& reduction_pos,
                unordered_map<SequenceType<StateType>, unsigned>& data)
{
    // modify the data map by erasing the previous keys and mapping the updated data with the new data count
    // after using the root node of aof to replace the leaf nodes, some of the data may collapse into one class

    vector<SequenceType<StateType>> keys; // store previous data
    unordered_map<SequenceType<StateType>, unsigned> new_data;

    for (auto d : data)
    {
        keys.push_back(d.first);
        SequenceType<StateType> new_key = d.first;

        vector<unsigned> positions = reduction_pos.at(d.first);
        sort(positions.begin(), positions.end()); 
        reverse(positions.begin(), positions.end());

        // replace aof from back to front
        for (auto p : positions)
        {
            new_key.erase(new_key.begin() + p, new_key.begin() + p + aof.or_children_.size());
            new_key.insert(new_key.begin() + p, aof.root_and_state_);
        }
        
        new_data[new_key] += d.second;
    }
}

template <class StateType, class AttributeType>
void Learn(AOG<StateType, AttributeType>& aog, unordered_map<SequenceType<StateType>, unsigned>& data,
           double alpha, int inner_loop_iters, int expand_limits)
{
    int outer_iters = 0;
    int outer_iter_limit = 50;
    int num_aofs = 0;
    unsigned next_state_id = aog.GetStates().size(); 
    Symbolic_State<StateType> root = aog.GetVertexContent(aog.GetRoot())->GetState();
        
    while (true)
    {
        cout << "Outer iteration " << outer_iters << "   Number of aofs selected: " << num_aofs << "\n";
        outer_iters++;

        vector<AOFStruct<StateType>> aof_set; // corresponding to F in the pseudocode
        vector<double> posterior_gains;
        posterior_gains.push_back(0); // use log gain instead of original gain
        
        double max_gain = 0;
        int max_gain_idx;
        int rd_grammar_total = 0;
        unsigned rd_total = 0;

        pair<unordered_map<SequenceType<StateType>, vector<unsigned>>, 
                 vector<unordered_map<Symbolic_State<StateType>, unsigned>>> reduction_for_max; // store the reduction for the aof with max gain

        for (int i = 0; i < inner_loop_iters; i++)
        {
            cout << "  Inner iteration " << i << "\n";

            Context_Matrix<StateType> cm;
            pair<SequenceType<StateType>, SequenceType<StateType>> init_extensions, extensions;

            // randomly generate a bigram
            // init_extensions stores possible leaf nodes at each position, which can be used later in the AddLeafs operation
            AOFStruct<StateType> aof = ChooseBigram(data, next_state_id, init_extensions);

            pair<unordered_map<SequenceType<StateType>, vector<unsigned>>, 
                 vector<unordered_map<Symbolic_State<StateType>, unsigned>>> reduction = GenerateCMandReduction(aof, data, cm, rd_total, rd_grammar_total, extensions);
                    
            max_gain = CalcPostGain(reduction, aof, cm, rd_total, rd_grammar_total, alpha);
            
            // expand the selected aof up to a maximum depth, use greedy search at each depth
            int max_depth = 1;
            while (max_depth)
            {
                max_gain_idx = -1;
                
                extensions.first.insert(extensions.first.end(), init_extensions.first.begin(), init_extensions.first.end());
                extensions.second.insert(extensions.second.end(), init_extensions.second.begin(), init_extensions.second.end());

                // obtain all possible expansions
                vector<AOFStruct<StateType>> expanded_aofs = ExpandAOFs(aof, extensions, next_state_id, expand_limits);
        
                for (int k = 0; k < expanded_aofs.size(); k++)
                {
                    AOFStruct<StateType> new_aof = expanded_aofs[k];
                    Context_Matrix<StateType> new_cm;

                    pair<unordered_map<SequenceType<StateType>, vector<unsigned>>, 
                            vector<unordered_map<Symbolic_State<StateType>, unsigned>>> new_reduction = GenerateCMandReduction(new_aof, data, new_cm, rd_total, rd_grammar_total, extensions);
                        
                    if (rd_total == 0)
                        continue;

                    double new_gain = CalcPostGain(new_reduction, new_aof, new_cm, rd_total, rd_grammar_total, alpha);

                    // if the new aof is better than the current optimal one, remember this new aof
                    if (new_gain > max_gain)
                    {
                        max_gain = new_gain;
                        max_gain_idx = k;
                        reduction = new_reduction;
                    }
                }

                max_depth--;

                if (max_gain_idx != -1)  
                    aof = expanded_aofs[max_gain_idx];
                else
                    break;

                extensions.first.clear();
                extensions.second.clear();
                GenerateCMandReduction(aof, data, cm, rd_total, rd_grammar_total, extensions);
            }
            
            // now we have the optimal aof generated from the current randomly selected bigram
            // check if it is the best among all randomly selected bigrams
            if (max_gain > posterior_gains.back())
            {
                aof_set.push_back(aof);
                posterior_gains.push_back(max_gain);
                reduction_for_max = reduction;
            }
            
        }

        if (aof_set.size() == 0)
            break;

        AOFStruct<StateType> selected_aof = aof_set.back();
        num_aofs++;
        
        // update data
        UpdateData(selected_aof, reduction_for_max.first, data);
        
        // add aof to aog
        vector<Symbolic_Rule<StateType>> new_rules = selected_aof.AOF2Rules();
        for (auto rule : new_rules)
            aog.AddRule(rule);

        // update weights within the aof in the aog
        for (int i = 0; i < selected_aof.or_children_.size(); i++)
        {
            unordered_map<VertexId, double> weights;
            VertexId or_child = aog.GetVertexIdByState(selected_aof.or_children_[i]);
            aog.SetIsNotOr(or_child, false); // set or for or_children
            
            for (auto cid : aog.ChildrenVertices(or_child))
            {
                // get the vertex representing the leafs
                Symbolic_State<StateType> cstate;
                if (aog.ChildrenVertices(cid).size() == 0)
                    cstate = aog.GetVertexContent(cid) -> GetState();
                else
                    cstate = aog.GetVertexContent(aog.ChildrenVertices(cid)[0]) -> GetState(); // because dummy nodes are added when new rules are added
                weights[cid] = reduction_for_max.second[i][cstate];
            }
            
            aog.SetOutEdgeWeights(or_child, weights);
        }

        // update aog by unlinking the leaf nodes to the data node and connecting the aof root node to the data node
        // first do the unlinking
        for(auto data_node : aog.ChildrenVertices(aog.GetRoot()))
        {
            vector<VertexId> children_ids = aog.ChildrenVertices(data_node);
            SequenceType<StateType> children;

            for(auto cid: children_ids)
                children.push_back(aog.GetVertexContent(cid) -> GetState()); // this represents the data instance
            
            aog.DeleteRule(Symbolic_Rule<StateType>(aog.GetVertexContent(aog.GetRoot()) -> GetState(), children));
        }

        // then do the reconnection
        for (auto d : data)
        {
            Symbolic_Rule<string> new_rule(aog.GetVertexContent(aog.GetRoot()) -> GetState(), d.first);
            aog.AddRule(new_rule);
        }

        // set weight for each data node in the aog to data count
        unordered_map<VertexId ,double> weights;
        for(auto data_node : aog.ChildrenVertices(aog.GetRoot()))
        {
            vector<VertexId> children_ids = aog.ChildrenVertices(data_node);
            SequenceType<StateType> children;

            for(auto cid : children_ids)
                children.push_back(aog.GetVertexContent(cid) -> GetState());

            weights[data_node] = data.at(children);
        }
        aog.SetOutEdgeWeights(aog.GetRoot(), weights);

        
        if (outer_iters >= outer_iter_limit || data.size() == 1)
            break;
    }
    return;
}