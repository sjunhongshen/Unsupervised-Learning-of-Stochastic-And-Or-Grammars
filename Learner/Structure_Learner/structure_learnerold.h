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
    //TODO

    cerr<<"\n^^^^^^^^^^^^^^^^^. current AOF\n";

    unordered_map<Symbolic_State<StateType>, double> freq_1;
    unordered_map<Symbolic_State<StateType>, double> freq_2;
    SequenceType<StateType> or1 = aof.leafs_[0];
    SequenceType<StateType> or2 = aof.leafs_[1];
    for (auto d : data)
    {
        for (int i = 0; i < d.first.size() - 1; i++)
        {
            // cerr<< "    Data size "<<d.first.size()<<endl;

            Symbolic_State<StateType> first = d.first[i];
            Symbolic_State<StateType> second = d.first[i + 1];
            if (find(or1.begin(), or1.end(), first) != or1.end() && find(or2.begin(), or2.end(), second) != or2.end())
            {
                freq_1[first] += d.second;
                freq_2[second] += d.second;
                // cerr<<"Find: "<< first.GetContent()<<" "<<second.GetContent()<<"   "<<d.second<<endl;
            }
        }
    }


    // cerr<<endl;
    for (auto l : freq_1)
    {
        cerr<<l.first.GetContent()<<" "<<l.second;
        // l.second = l.second / appearance_count;
        cerr<<endl;
    }
    for (auto l : freq_2)
    {
        cerr<<l.first.GetContent()<<" "<<l.second;
        // l.second = l.second / appearance_count;
        cerr<<endl;
    }

    vector<SequenceType<StateType>> keys;
    unordered_map<SequenceType<StateType>, unsigned> new_data;

    for (auto d : data)
    {
        keys.push_back(d.first);

        vector<unsigned> positions = reduction_pos.at(d.first);

        sort(positions.begin(), positions.end()); 
        reverse(positions.begin(), positions.end());

        SequenceType<StateType> new_key = d.first;
        // cerr<<"\n\n        UPDATED data  BEFORE. ";
        // for (auto l : new_key)
        // {
        //     cerr<<l.GetContent()<<" ";
        // }

        for (auto p : positions)
        {
            new_key.erase(new_key.begin() + p, new_key.begin() + p + aof.or_children_.size());
            new_key.insert(new_key.begin() + p, aof.root_and_state_);//?????
        }

        // cerr<<"\n        UPDATED data  AFTER ";
        // for (auto l : new_key)
        // {
        //     cerr<<l.GetId()<<" ";
        // }
        
        new_data[new_key] += d.second;
    }
    double sum = 0;
    for (int i = 0; i < keys.size(); i++)
    {
        sum += data[keys[i]];
        // cerr<<" old count "<< data[keys[i]]<<" "<< sum<<endl;
        data.erase(keys[i]);

    }
    sum = 0;
    for (auto d : new_data)
    {
        data[d.first] = d.second;
        sum += d.second;
        // cerr<<" new count "<< d.second<<" "<< sum<<endl;
    }
    

}

template <class StateType, class AttributeType>
void Learn(AOG<StateType, AttributeType>& aog, unordered_map<SequenceType<StateType>, unsigned>& data,
           double alpha, int inner_loop_iters, int expand_limits)
{
    int outer_iters = 0;
    int outer_iter_limit = 50;
    unsigned next_state_id = aog.GetStates().size();
    unordered_map<VertexId, unordered_map<VertexId, double>> ors;
    // cerr<<"next state id "<< next_state_id<<endl;

    // for (auto s: aog.GetStates())
    // {
    //     cerr<<"state "<< s.GetContent()<<endl;
    // }

    // for(auto data_node : aog.ChildrenVertices(aog.GetRoot()))
    // {
    //     cerr<<data_node<<endl;
    //     vector<VertexId> children_ids = aog.ChildrenVertices(data_node);

    //     for(auto cid: children_ids)
    //     {
    //         aog.AddRule(Symbolic_Rule<StateType>(aog.GetVertexContent(data_node)->GetState(), {aog.GetVertexContent(cid)->GetState()}));
    //     }
    // }
    
    Symbolic_State<StateType> root = aog.GetVertexContent(aog.GetRoot())->GetState();
    int num_aofs = 0;
    double prev_gain = -INFINITY;

    // cerr<<"ROOT "<< root.GetContent()<<endl;
    
    while(true)
    {
        outer_iters++;
        //TODO
        vector<AOFStruct<StateType>> aof_set; //F
        vector<double> posterior_gains;
        posterior_gains.push_back(0);
        
        
        double max_gain;
        int max_gain_idx;
        int rd_grammar_total = 0;
        unsigned rd_total = 0;

        pair<unordered_map<SequenceType<StateType>, vector<unsigned>>, 
                 vector<unordered_map<Symbolic_State<StateType>, unsigned>>> reduction_for_max;

        for (int i = 0; i < inner_loop_iters; i++)
        {
            

            Context_Matrix<StateType> cm;
            pair<SequenceType<StateType>, SequenceType<StateType>> init_extensions, extensions;

            AOFStruct<StateType> aof = ChooseBigram(data, next_state_id, init_extensions);

            pair<unordered_map<SequenceType<StateType>, vector<unsigned>>, 
                 vector<unordered_map<Symbolic_State<StateType>, unsigned>>> reduction = GenerateCMandReduction(aof, data, cm, rd_total, rd_grammar_total, extensions);
                    
            if (rd_total == 0)
                continue;
            max_gain = CalcPostGain(reduction, aof, cm, rd_total, rd_grammar_total, alpha);
            

            cerr <<"******************* EXPAND RANDOM AOF outer round "<<outer_iters<< "   inner round. "<<i<<"  orginal max gain "<<max_gain<<"\n";

            int max_depth = 3;
            while (max_depth)
            {
                max_gain_idx = -1;
                
                extensions.first.insert(extensions.first.end(), init_extensions.first.begin(), init_extensions.first.end());
                extensions.second.insert(extensions.second.end(), init_extensions.second.begin(), init_extensions.second.end());

                vector<AOFStruct<StateType>> expanded_aofs = ExpandAOFs(aof, extensions, next_state_id, expand_limits);

                extensions.first.clear();
                extensions.second.clear();
        
                for (int k = 0; k < expanded_aofs.size(); k++)
                {
                    // cerr <<"      ******************* EXPAND "<<k<<"  of total "<< expanded_aofs.size()<<endl;

                    AOFStruct<StateType> new_aof = expanded_aofs[k];
                    Context_Matrix<StateType> new_cm;

                    pair<unordered_map<SequenceType<StateType>, vector<unsigned>>, 
                            vector<unordered_map<Symbolic_State<StateType>, unsigned>>> new_reduction = GenerateCMandReduction(new_aof, data, new_cm, rd_total, rd_grammar_total, extensions);
                        
                    if (rd_total == 0)
                        continue;
                    double new_gain = CalcPostGain(new_reduction, new_aof, new_cm, rd_total, rd_grammar_total, alpha);

                    if (new_gain > max_gain)
                    {
                        max_gain = new_gain;
                        max_gain_idx = k;
                        reduction = new_reduction;
                    }
                    // cerr <<"      ******************* END EXPAND   gain: "<<new_gain<<" prev max gain "<<max_gain<<" idx "<<max_gain_idx<<endl;
                }

                max_depth--;

                if (max_gain_idx != -1)      
                    aof = expanded_aofs[max_gain_idx];
                else
                    break;

            }
            
            cerr <<"******************* END EXPAND RANDOM AOF  \n";
            if (max_gain > posterior_gains.back())
            {
                aof_set.push_back(aof);
                posterior_gains.push_back(max_gain);
                reduction_for_max = reduction;
                cerr <<"******************* Greater than before   "<<max_gain<<"   total num of aof "<<num_aofs<<"   \n";
            }
            

        }

        if (aof_set.size() == 0)
            break;

        AOFStruct<StateType> selected_aof = aof_set.back();
        prev_gain = posterior_gains.back();
        cerr << "!! BEST "<<prev_gain<<endl;
        // Context_Matrix<StateType> cm;
        // pair<SequenceType<StateType>, SequenceType<StateType>> extensions;//?


        // pair<unordered_map<SequenceType<StateType>, vector<unsigned>>, 
        //     vector<unordered_map<Symbolic_State<StateType>, unsigned>>> reduction = GenerateCMandReduction(selected_aof, data, cm, rd_total, rd_grammar_total, extensions);

        UpdateData(selected_aof, reduction_for_max.first, data);
        num_aofs++;
        vector<Symbolic_Rule<StateType>> new_rules = selected_aof.AOF2Rules();

        cerr<<":::::::::::::::: data size "<<data.size()<<"  grammar size. "<<aog.NumOfRules()<< endl;
        int idx = 0;
        for (auto rule : new_rules)
        {
            aog.AddRule(rule);
            idx++;
            // break;
        }
        cerr<<":::::::::::::::: data size "<<data.size()<<"  grammar size. "<<aog.NumOfRules()<< endl;

        for (int i = 0; i < selected_aof.or_children_.size(); i++)
        {
            unordered_map<VertexId, double> weights;
            VertexId or_child = aog.GetVertexIdByState(selected_aof.or_children_[i]);
            // cerr<<"@@@@@@@@@@@@@@  orid  "<<or_child<<". "<<aog.GetVertexContent(or_child)->GetState().GetContent()<<". "<< selected_aof.or_children_[i].GetId()<<endl;
            aog.SetIsNotOr(or_child, false);
            
            // VertexId or_child = aog.GetVertexIdByState(selected_aof.or_children_[i]);
            // cerr<<"@@@@@@@@@@@@@@  orid  "<<or_child<<". "<<aog.GetVertexContent(or_child)->GetState().GetContent()<<endl;

            for (auto cid : aog.ChildrenVertices(or_child))
            {
                // cerr<<"@@@@@@@@@@@@@@  cid  "<<cid<<". "<<aog.GetVertexContent(cid)->GetState().GetId()<<endl;
                Symbolic_State<StateType> cstate;
                if (aog.ChildrenVertices(cid).size() == 0)
                    cstate = aog.GetVertexContent(cid)->GetState();
                else
                    cstate = aog.GetVertexContent(aog.ChildrenVertices(cid)[0])->GetState();
                weights[cid] = reduction_for_max.second[i][cstate];
            }
            // for (int j = 0; j < selected_aof.leafs_[i].size(); j++)
            // {
            //     VertexId dummy = aog.ParentsVertices(aog.GetVertexIdByState(selected_aof.leafs_[i][j]))[0];
            //     cerr<<"dummy     "<<dummy<<endl;
            //     weights[dummy] = reduction_for_max.second[i][selected_aof.leafs_[i][j]];
            //     cerr<<"&&&&&&&&&    set weight for "<<selected_aof.leafs_[i][j].GetContent()<<"  id "<<aog.GetVertexIdByState(selected_aof.leafs_[i][j])<<"  "<<reduction_for_max.second[i][selected_aof.leafs_[i][j]]<<" "<< i <<" "<<j<<endl;
            // }
            ors[or_child] = weights;
            aog.SetOutEdgeWeights(or_child, weights);
        }

        // unordered_map<VertexId ,double> weights = aog.GetOutEdgeWeights(aog.GetRoot(), false);
        // for(auto it = weights.begin(); it != weights.end(); it++)
        // {
            
        //     cerr << "weights "<<it->first<<" "<<weights[it->first]<<endl;
        // }
        // for (auto d: aog.GetRulesAsSource(aog.GetVertexContent(aog.GetRoot())->GetState()))
        //     for (auto r:d.GetResults())
        //         cerr<<"WTF. "<<r.GetContent()<<endl;
        for(auto data_node : aog.ChildrenVertices(aog.GetRoot()))
        {
            // cerr<<data_node<<endl;
            vector<VertexId> children_ids = aog.ChildrenVertices(data_node);
            SequenceType<StateType> children;

            for(auto cid: children_ids)
            {

                children.push_back(aog.GetVertexContent(cid)->GetState());
                // cerr<<"&&&&&&&&&&.    children "<<cid<<" "<<aog.GetVertexContent(cid)->GetState().GetContent()<<endl;
            }
            
            // weights[it->first] = rules_count.at(children);
            if(!aog.DeleteRule(Symbolic_Rule<StateType>(aog.GetVertexContent(aog.GetRoot())->GetState(), children)))
                cerr<<"WTF"<<endl;
            // vector<unsigned> pos = reduction_for_max.first[children];
        //     for (auto p : pos)
        //     {
        //         for (int i = 0; i < selected_aof.or_children_.size(); i++)
        //         {
        //             if(!aog.DeleteRule(Symbolic_Rule<StateType>(aog.GetVertexContent(data_node)->GetState(), {children[p + i]})))
        //             {
        //                 // if(!aog.DeleteRule(Symbolic_Rule<StateType>(aog.GetVertexContent(aog.GetRoot())->GetState(), {aog.GetVertexContent(data_node)->GetState()})))
        //                 //     cerr<<"WTF"<<endl;
        //                 cerr<<"--- "<< p<<" "<<children[p + i].GetContent()<<endl;
        //             }
        //         }
        // //         
        // //         cerr<<"root and state. "<<aog.GetVertexIdByState(selected_aof.root_and_state_)<<endl;
        //     }
        //     aog.AddRule(Symbolic_Rule<StateType>(aog.GetVertexContent(data_node)->GetState(), {selected_aof.root_and_state_}));
            // cerr<<"root and state. "<<aog.GetVertexIdByState(selected_aof.root_and_state_)<<endl;
        }
        for (auto d: data)
        {
            // rules_count[dataset[i]]++;
            // cerr<<dataset[i][0].GetContent()<<endl;
            Symbolic_Rule<string> new_rule(aog.GetVertexContent(aog.GetRoot())->GetState(), d.first);
            aog.AddRule(new_rule);
        }
        unordered_map<VertexId ,double> weights;
        for(auto data_node : aog.ChildrenVertices(aog.GetRoot()))
        {
            // cerr<<data_node<<endl;
            vector<VertexId> children_ids = aog.ChildrenVertices(data_node);
            SequenceType<StateType> children;

            for(auto cid: children_ids)
            {

                children.push_back(aog.GetVertexContent(cid)->GetState());
                // cerr<<"&&&&&&&&&&.    children "<<cid<<" "<<aog.GetVertexContent(cid)->GetState().GetContent()<<endl;
            }
            
            weights[data_node] = data.at(children);

        }
        aog.SetOutEdgeWeights(aog.GetRoot(), weights);

        
        // for(auto it = weights.begin(); it != weights.end(); it++)
        // {
        //     vector<VertexId> children_id = aog.ChildrenVertices(it->first);
        //     SequenceType<string> children;
        //     for(auto cid: children_id)
        //     {
        //         children.push_back(aog.GetVertexContent(cid)->GetState());
        //         cerr<<aog.GetVertexContent(cid)->GetState().GetContent();
        //     }
        //     weights[it->first] = data.at(children);
        //     cerr << "weights "<<it->first<<" "<<weights[it->first]<<endl;
        // }
        // aog.SetOutEdgeWeights(aog.GetRoot(), weights);
        cerr<<":::::::::::::::: data size "<<data.size()<<"  grammar size. "<<aog.NumOfRules()<< endl;

        // for (auto o :ors)
        // {
        //     aog.SetIsNotOr(o.first, false);
        //     aog.SetOutEdgeWeights(o.first, o.second);
        // }
        if(outer_iters >= outer_iter_limit || data.size() == 1)
            break;
    }
    return;
}