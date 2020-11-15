#include <vector>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <math.h>
#include <algorithm>
#include <numeric>

#include "AOG.h"
#include "Context_Matrix.h"

using namespace std;
using namespace AOG_LIB;

std::mt19937 gen;

template<class StateType>
struct AOFStruct
{
    Symbolic_State<StateType> root_and_state_;
    SequenceType<StateType> or_children_;
    vector<SequenceType<StateType>> leafs_;
    vector<unordered_map<Symbolic_State<StateType>, double> > weights_;

    // return the rules defining this And-Or fragment
    vector<Symbolic_Rule<StateType> > AOF2Rules() const
    {
        assert(or_children_.size() == leafs_.size());
        vector<Symbolic_Rule<StateType> > rules;
        rules.push_back(Symbolic_Rule<string>(this->root_and_state_, this->or_children_));
        for(int i = 0; i < this->or_children_.size(); i++)
        {
            for(int j = 0; j < this->leafs_[i].size();j++)
                rules.push_back(Symbolic_Rule<string>(this->or_children_[i], {this->leafs_[i][j]}));
        }
        
        return rules;
    }

    // return all possible configurations of this And-Or fragment
    // configurations are all combination of the or-node leafs
    unordered_set<SequenceType<StateType>> GetConfigs()
    {
        unordered_set<SequenceType<StateType>> all_configs;
        stack<SequenceType<StateType>> result_stack;
        if(this->leafs_.size() == 1)
        {
            for(int i = 0; i < this->leafs_[0].size(); i++)
                all_configs.insert({this->leafs_[0][i]});
        }
        else
        {
            for(int i = 0; i < this->leafs_[0].size(); i++)
                result_stack.push({this->leafs_[0][i]});
            while(!result_stack.empty())
            {
                SequenceType<StateType> temp = result_stack.top();
                result_stack.pop();
                for(int i = 0; i < this->leafs_[temp.size()].size(); i++)
                {
                    temp.push_back(this->leafs_[temp.size()][i]);
                    if(temp.size() == this->or_children_.size())
                        all_configs.insert(temp);
                    else
                        result_stack.push(temp);
                    temp.pop_back();
                }
            }
        }
        return all_configs;
    }
};


/**
 * randomly generate an And-Or fragment out of current data
 * 
 * @param data: mapping data to the number of its occurance, the occurance is also the weight from the root node
 * @param state_idx: the index of the new non-terminal states created in this function
 * 
 * @returns an And-Or fragment with two-or nodes and two leafs for each or-node,
 *          empty if there cannot be such an AOF generated from current dataset
 */ 
template<class StateType>
AOFStruct<StateType> ChooseBigram(unordered_map<SequenceType<StateType>, unsigned>& data, unsigned& state_idx, pair<SequenceType<StateType>, SequenceType<StateType>>& extensions)
{
    AOFStruct<StateType> aof;
    Symbolic_State<StateType> and_root(state_idx++);
    Symbolic_State<StateType> or_child_1(state_idx++);
    Symbolic_State<StateType> or_child_2(state_idx++);
    aof.root_and_state_ = and_root;
    aof.or_children_.push_back(or_child_1);
    aof.or_children_.push_back(or_child_2);

    // randomly generate two numbers to represent the data indices
    std::uniform_int_distribution<> dis1(0, data.size() - 1);
    int data_idx1 = dis1(gen);
    int data_idx2 = data_idx1;
    while (data_idx2 == data_idx1)
        data_idx2 = dis1(gen);

    int idx = 0;
    SequenceType<StateType> or1, or2;

    for (auto d : data)
    {
        if(idx == data_idx1 || idx == data_idx2)
        {
            if (d.first.size() == 1)
            {
                or1.push_back(d.first[0]); // if the chosen data only contain one state
                idx++;
                continue;
            }

            // randomly generate a position index
            std::uniform_int_distribution<> dis2(0, d.first.size() - 2);
            int node_idx = dis2(gen);
            
            // we don't want the two or children to contain the same state because it is rare for two same words to appear consecutively in a sentence
            // though in reality we may encounter such sentences, but it's rare
            int max_count = 10;
            while (max_count > 0 && (or1.size() != 0 && d.first[node_idx + 1] == or1.back() || or2.size() != 0 && d.first[node_idx] == or2.back())) // (d.first[node_idx].GetId()!= -1 || d.first[node_idx + 1].GetId()!= -1)
            {
                node_idx = dis2(gen);
                max_count--;
            }

            if (or1.size() == 0 || d.first[node_idx] != or1.back())
                or1.push_back(d.first[node_idx]);
            if (or2.size() == 0 || d.first[node_idx + 1] != or2.back())
                or2.push_back(d.first[node_idx + 1]);    
        }
        idx++;
    }

    aof.leafs_.push_back(or1);
    aof.leafs_.push_back(or2);

    // count number of appearance of the aof in the entire dataset
    // meanwhile look for possible candidate leaf nodes to be added in AddLeafs operation
    unordered_map<Symbolic_State<StateType>, double> freq_1;
    unordered_map<Symbolic_State<StateType>, double> freq_2;
    for (auto d : data)
    {
        for (int i = 0; i < d.first.size() - 1; i++)
        {
            Symbolic_State<StateType> first = d.first[i];
            Symbolic_State<StateType> second = d.first[i + 1];

            // config found
            if (find(or1.begin(), or1.end(), first) != or1.end() && find(or2.begin(), or2.end(), second) != or2.end())
            {
                freq_1[first] += d.second;
                freq_2[second] += d.second;
            }
            else if (find(or1.begin(), or1.end(), first) != or1.end())
            {
                extensions.second.push_back(second);
            }
            else if (find(or2.begin(), or2.end(), second) != or2.end())
            {
                extensions.first.push_back(first);
            }
        }
    }

    aof.weights_.push_back(freq_1);
    aof.weights_.push_back(freq_2);

    return aof;
}


/**
 * recursively search for appearance of configs in a piece of data
 *
 * @param configs: possible configs of an aof, target to find
 * @param current_config: existing configs in the data, we add to this if we find a config
 * @param context: context of the data, we add to this while searching the data
 * @param length: number of or children in the aof
 * @param pos: we add the starting index to this vector if we find a config
 *
 */ 
template<class StateType>
void search_data(unordered_set<SequenceType<StateType>>& configs, vector<SequenceType<StateType> >& current_config, ContextType<StateType>& context, SequenceType<StateType>& data, int length, unsigned start_pos, vector<unsigned>& pos)

{
    int find = 0;
    if (data.size() < length) // data size smaller than the aof size
    {
        context.push_back(data);
        return;
    }

    for (int i = 0; i <= data.size() - length; i++)
    {
        for (auto config : configs)   
        {
            SequenceType<StateType> partial_data(data.begin() + i, data.begin() + i + length);

            if (partial_data == config)
            {
                current_config.push_back(config);
                pos.push_back(start_pos + i);

                // context before the config
                if (i != 0)
                {
                    SequenceType<StateType> before(data.begin(), data.begin() + i);
                    context.push_back(before);
                }

                // search the rest of the data
                SequenceType<StateType> after(data.begin() + i + length, data.end());
                search_data(configs, current_config, context, after, length, start_pos + i + length, pos);
                find = 1;
                break;
            }
        }
        if (find)
            break;
    }

    if (!find)
        context.push_back(data);
}

/**
 * for a given And-Or fragment, generate its context matrix and reductions
 * weights of the AOF is also updated in this function
 * 
 * @param aof: the query and-or fragment
 * @param data: mapping data to the number of its occurance, the occurance is also the weight from the root node
 * @param cm: the context matrix, updated in the function pass by reference
 * @param rd_total: number of reductions introduced by this aof,
 *                  need to count for each data (consider data occurance)
 * @param rd_grammar_total: number of reductions in the grammar rules,
 *                          NO need to count for each data (NO consider data occurance), used for prior gain
 * @param extensions: consecutive states (before and after) show up next to the aof configurations
 *                    used for the four operations to generate new aofs
 * 
 * @returns reduction of this aof, a pair (first: data->replacement positions, used for data update;
 *                                         second: state->occurance, used for or-node weight calculation)
 */ 
template<class StateType>
pair<unordered_map<SequenceType<StateType>, vector<unsigned>>, vector<unordered_map<Symbolic_State<StateType>, unsigned>>>
GenerateCMandReduction(AOFStruct<StateType>& aof,
                       const unordered_map<SequenceType<StateType>, unsigned>& data,
                       Context_Matrix<StateType>& cm, unsigned& rd_total, int& rd_grammar_total,
                       pair<SequenceType<StateType>, SequenceType<StateType>>& extensions)
{
    assert(aof.or_children_.size() > 0);
    rd_total = 0;
    rd_grammar_total = 0;

    pair<unordered_map<SequenceType<StateType>, vector<unsigned>>,
         vector<unordered_map<Symbolic_State<StateType>, unsigned>>> reduction;
    
    // initialize state -> occurance map for each or child
    for (int i = 0; i < aof.or_children_.size(); i++)
    {
        unordered_map<Symbolic_State<StateType>, unsigned> count;
        for (int j = 0; j < aof.leafs_[i].size(); j++)
            count[aof.leafs_[i][j]] = 0;
        reduction.second.push_back(count); 
    }

    unordered_set<SequenceType<StateType>> configs = aof.GetConfigs();

    for (auto d : data)
    {
        SequenceType<StateType> data_copy = d.first;        
        ContextType<StateType> context;
        vector<SequenceType<StateType> > current_config;
        vector<unsigned> pos;

        // search for possible configs
        search_data(configs, current_config, context, data_copy, aof.or_children_.size(), 0, pos);
        reduction.first[d.first] = pos;

        for (auto p : pos)
        {
            if (p >= 1)
                extensions.first.push_back(d.first[p - 1]); // possible states for adding or-nodes before the current aof

            if (p + aof.or_children_.size() < d.first.size())
                extensions.second.push_back(d.first[p + aof.or_children_.size()]); // possible states for adding or-nodes after the current aof
        }

        if (current_config.size() > 0)
        {
            auto elements_by_context = boost::multi_index::get<0>(cm).equal_range(context);
            int add_new = 1;

            // see if we already have an entry in the context matrix that is exactly the same
            for (auto &element = elements_by_context.first; element != elements_by_context.second; element++)
            {
                if (element -> configuration_ == current_config)
                {
                    add_new = 0;
                    for (int j = 0; j < d.second; j++)
                        boost::multi_index::get<0>(cm).modify(element, [](AOG_LIB::CM_Entry<StateType> &element) {element.count_++; }); // increament count
                }
            }

            // no existing entry, add a new one
            if (add_new)
                cm.insert({context, current_config, d.second});

            // update reduction map, rd_total and rd_grammar_total
            for (auto config : current_config)
            {
                for (int k = 0; k < aof.or_children_.size(); k++)
                {
                    reduction.second[k][config[k]] += d.second; 
                }
                rd_total += d.second;
                rd_grammar_total += 1;
            }
        }
    }

    return reduction;
}


// Calculate the posterior gain of an And-Or fragment
template<class StateType>
double CalcPostGain(const pair<unordered_map<SequenceType<StateType>, vector<unsigned>>,
                               vector<unordered_map<Symbolic_State<StateType>, unsigned>>>& reduction,
                    const AOFStruct<StateType>& aof, const Context_Matrix<StateType>& cm,
                    const unsigned rd_total, const int rd_grammar_total, double alpha)
{
    // all probability gains are calculated in log scale

    double post_gain = 0, likelihood_gain = 0;
    double children_num = 0;
    
    // first calculate likelihood gain
    
    // first part of likelihood gain formula
    for (int i = 0; i < aof.or_children_.size(); i++)
    {
        for (int j = 0; j < aof.leafs_[i].size(); j++)
        {
            Symbolic_State<StateType> l(aof.leafs_[i][j]);
            likelihood_gain += reduction.second[i].at(l) * log((((double)reduction.second[i].at(l)) + 1e-7) / rd_total);

            if (reduction.second[i].at(l) != 0)
                children_num++;
        }
    }

    // index context matrix by context    
    auto& elements_by_contexts = boost::multi_index::get<0>(cm);
    vector<vector<SequenceType<StateType> >> counted_contexts;

    for (auto same_context = elements_by_contexts.begin(); same_context != elements_by_contexts.end(); same_context++)
    {
        // if the context has already been calculated, continue
        if (find(counted_contexts.begin(), counted_contexts.end(), same_context -> context_) != counted_contexts.end())
            continue;
        
        auto element = elements_by_contexts.equal_range(same_context -> context_); // all entries with the same context
        
        double sum = 0;
        for (auto e = element.first; e != element.second; e++)
        {
            sum += e -> count_;
            likelihood_gain -= e -> count_ * log(e -> count_ + 1e-7); // denominator of the second part of the formula
        }

        likelihood_gain += sum * log(sum + 1e-7); // numerator of the second part
        counted_contexts.push_back(same_context -> context_);   
    }

    // now calculate the prior gain
    double prior_gain = -alpha * (aof.or_children_.size() + children_num - rd_grammar_total * (aof.or_children_.size() - 1));
    
    post_gain = likelihood_gain + prior_gain; // in log scale

    return post_gain;
}

// Four types of Aod-Or Fragment editing operations
template<class StateType>
vector<AOFStruct<StateType>> DeleteLeaf(const AOFStruct<StateType>& aof)
{
    vector<AOFStruct<StateType>> new_aofs;
    for(int i = 0; i < aof.or_children_.size(); i++)
    {
        if(aof.leafs_[i].size() == 1)
            continue;
        for(int j = 0; j < aof.leafs_[i].size(); j++)
        {
            AOFStruct<StateType> temp = aof;
            temp.leafs_[i].erase(temp.leafs_[i].begin() + j);
            new_aofs.push_back(temp);
        }
    }
    return new_aofs;
}

template<class StateType>
vector<AOFStruct<StateType>> DeleteOr(const AOFStruct<StateType>& aof)
{
    vector<AOFStruct<StateType>> new_aofs;
    if(aof.or_children_.size() == 1)
        return new_aofs;
    for(int i = 0; i < aof.or_children_.size(); i++)
    {
        AOFStruct<StateType> temp = aof;
        temp.or_children_.erase(temp.or_children_.begin() + i);
        temp.leafs_.erase(temp.leafs_.begin() + i);
        new_aofs.push_back(temp);
    }

    return new_aofs;
}

// only add leafs at first and second or-nodes
template<class StateType>
vector<AOFStruct<StateType>> AddLeafs(const AOFStruct<StateType>& aof,
                                      const SequenceType<StateType>& states_ahead,
                                      const SequenceType<StateType>& states_behind)
{
    vector<AOFStruct<StateType>> new_aofs;
    for(auto state: states_ahead)
    {
        bool repetitive = false;
        for(auto leaf: aof.leafs_[0])
            if(state == leaf)
            {
                repetitive = true;
                break;
            }
        if(repetitive)
            continue;
        AOFStruct<StateType> temp = aof;
        temp.leafs_[0].push_back(state);
        new_aofs.push_back(temp);
    }
    for(auto state: states_behind)
    {
        bool repetitive = false;
        for(auto leaf: aof.leafs_.back())
            if(state == leaf)
            {
                repetitive = true;
                break;
            }
        if(repetitive)
            continue;
        AOFStruct<StateType> temp = aof;
        temp.leafs_.back().push_back(state);
        new_aofs.push_back(temp);
    }
    return new_aofs;
}

// only add an or-node at the beginning or end
template<class StateType>
vector<AOFStruct<StateType>> AddOrs(const AOFStruct<StateType>& aof,
                                    const SequenceType<StateType>& states_ahead,
                                    const SequenceType<StateType>& states_behind,
                                    unsigned& state_idx, int expand_limits)
{
    vector<AOFStruct<StateType>> new_aofs;
    vector<int> states_ahead_indices(states_ahead.size());
    iota(states_ahead_indices.begin(), states_ahead_indices.end(), 0);
    vector<int> states_behind_indices(states_behind.size());
    iota(states_behind_indices.begin(), states_behind_indices.end(), 0);
    if(states_ahead.size() > expand_limits)
    {
        random_shuffle(states_ahead_indices.begin(), states_ahead_indices.end());
        states_ahead_indices.resize(int(sqrt(2 * expand_limits)));
    }
    if(states_behind.size() > expand_limits)
    {
        random_shuffle(states_behind_indices.begin(), states_behind_indices.end());
        states_behind_indices.resize(int(sqrt(2 * expand_limits)));
    }
    for(int i = 0; i < states_behind_indices.size(); i++)
    {
        for(int j = i + 1; j < states_behind_indices.size(); j++)
        {
            AOFStruct<StateType> temp = aof;
            Symbolic_State<StateType> new_state(state_idx++);
            temp.or_children_.push_back(new_state);
            temp.leafs_.push_back({states_behind[states_behind_indices[i]], states_behind[states_behind_indices[j]]});
            new_aofs.push_back(temp);
        }
    }
    for(int i = 0; i < states_ahead_indices.size(); i++)
    {
        for(int j = i + 1; j < states_ahead_indices.size(); j++)
        {
            AOFStruct<StateType> temp = aof;
            Symbolic_State<StateType> new_state(state_idx++);
            temp.or_children_.insert(temp.or_children_.begin(), new_state);
            temp.leafs_.insert(temp.leafs_.begin(), {states_ahead[states_ahead_indices[i]], states_ahead[states_ahead_indices[j]]});
            new_aofs.push_back(temp);
        }
    }

    return new_aofs;
}

template<class StateType>
vector<AOFStruct<StateType>> ExpandAOFs(const AOFStruct<StateType>& aof,
                                        const pair<SequenceType<StateType>, SequenceType<StateType> >& extension,
                                        unsigned& state_idx, int expand_limits)
{
    vector<AOFStruct<StateType>> new_aofs;
    vector<AOFStruct<StateType>> temp = DeleteLeaf(aof);
    new_aofs.insert(new_aofs.end(), temp.begin(), temp.end());
    temp = DeleteOr(aof);
    new_aofs.insert(new_aofs.end(), temp.begin(), temp.end());
    int position = new_aofs.size();
    //remove duplicate states
    unordered_set<Symbolic_State<StateType>> states_ahead_set;
    for(auto s: extension.first)
        states_ahead_set.insert(s);
    unordered_set<Symbolic_State<StateType>> states_behind_set;
    for(auto s: extension.second)
        states_behind_set.insert(s);
    SequenceType<StateType> states_ahead;
    for(auto it = states_ahead_set.begin(); it != states_ahead_set.end(); it++)
        states_ahead.push_back(*it);
    SequenceType<StateType> states_behind;
    for(auto it = states_behind_set.begin(); it != states_behind_set.end(); it++)
        states_behind.push_back(*it);

    temp = AddLeafs(aof, states_ahead, states_behind);
    new_aofs.insert(new_aofs.end(), temp.begin(), temp.end());
    temp = AddOrs(aof, states_ahead, states_behind, state_idx, expand_limits);
    new_aofs.insert(new_aofs.end(), temp.begin(), temp.end());
    if(new_aofs.size() > expand_limits)
    {
        random_shuffle(new_aofs.begin() + position, new_aofs.end());
        new_aofs.resize(expand_limits);
    }
    return new_aofs;
}