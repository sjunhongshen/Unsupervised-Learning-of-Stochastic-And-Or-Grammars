//
// Created by yuanluyao on 10/23/17.
//

#ifndef AOG_LIB_AOG_H
#define AOG_LIB_AOG_H

#include <memory>
#include <random>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <iostream>
#include <stack>
#include <utility>
#include <algorithm>
#include <chrono>
#include <math.h>
typedef std::chrono::high_resolution_clock Clock;

#include "../Core/Graph.hpp"
#include "AOG_Vertex.h"
#include "AOG_Edge.h"
#include "Symbolic_Rule.h"

template <class StateType, class AttributeType>
using AOG_Graph = AOG_LIB::Graph<AOG_LIB::AOG_Vertex<StateType, AttributeType>, AOG_LIB::AOG_Edge>;

namespace AOG_LIB
{
    /**
     * This class defines an AOG graph with basic operations
     * An AOG graph consists of And-nodes and Or-nodes with 
     * their own semantic meanings and edges with different weights.
     * Each vertex in the AOG graph has its own unique id.
     */
    template<class StateType, class AttributeType>
    class AOG : public AOG_Graph<StateType, AttributeType>
    {

        /**
         * all the rules stored in a hash table for easy look up
         */
        std::unordered_set<Symbolic_Rule<StateType> > all_rules_;
        /**
         * all leaf states stored in a hash table for easy look up
         */
        std::unordered_set<Symbolic_State<StateType> > all_leaf_states_;
        /**
         * a map from state to graph vertex id
         */
        std::unordered_map<Symbolic_State<StateType>, VertexId> state_to_vertex_;
        /**
         * the vertex id of the root node
         */
        VertexId root_id_;
        /**
         * an auxiliary boolean variable for error checking when adding a vertex
         */
        bool has_root_;

        /**
         * @brief: map vertexId pairs to their potential function pointers
         * 
         */
        std::unordered_map<std::pair<VertexId, VertexId>,
                           std::function<double(const std::vector<std::vector<AttributeType> >&, void*)>,
                           pair_hash> potential_funcs_;

        /**
         * @brief: map vertexId pairs to their gradient function
         * of their potential function pointers
         * 
         */
        std::unordered_map<std::pair<VertexId, VertexId>,
                           std::function<void(const std::vector<std::vector<AttributeType> >&, void*, std::vector<double>&)>,
                           pair_hash> potential_grad_funcs_;
        
        /**
         * @brief: map vertexId pairs to their potential function parameters
         * 
         */
        std::unordered_map<std::pair<VertexId, VertexId>, void*, pair_hash> potential_params_;

    private:
        /**
         * This function adds an edge between two vertices.      
         * Cannot add an edge that starts and ends at same vertex
         * AOG object should call AddRule instead of AddEdge
         * 
         * @param source: the vertex id of the source node
         * @param target: the vertex id of the target node
         * @param aog_edge: the edge to be added from the source node to the target node
         * @param multi_edge: true if the AOG graph allows multiple edges between two vertices
         * 
         * @returns true if the edge is added successfully
         */        
        virtual bool AddEdge(const VertexId, const VertexId,
                             const std::shared_ptr<AOG_Edge>,
                             bool = true, int = -1);

        /**
         * Disable DeleteEdge function for AOG class, use DeleteRule instead
         */
        using AOG_Graph<StateType, AttributeType>::DeleteEdge;

    public:
        /**
         * default constructor
         */
        AOG();
    
        /**
         * construct from a set of rules. Currently it assigns the weight of 1.0 to all edges
         */
        AOG(const std::vector<Symbolic_Rule<StateType> >&);

        /**
         * AOG copy constructor
         */
        AOG(const AOG& other);

        /**
         * AOG constructor by reading in a file
         */
        AOG(const std::string&);

        /**
         * construct from a set of leaf states
         */
        explicit AOG(const std::vector<Symbolic_State<StateType> > &);
        /**
         * @returns nRules the number of rules the AOG contains
         */
        unsigned NumOfRules() const { return unsigned(this->all_rules_.size()); }
        /**
         * @returns nSymbolic the number of symbolic states in the AOG
         */
        unsigned NumOfStates() const { return unsigned(this->state_to_vertex_.size()); }
        /**
         * @returns nLeaftStates the number of leaf states in the AOG
         */
        unsigned NumOfLeafStates() const { return unsigned(this->all_leaf_states_.size()); }
        /**
         * @returns nVertices the number of vertices in the AOG
         */
        virtual unsigned NumberOfVertices() const {return state_to_vertex_.size();}
        /**
         * @returns is_root whether the root of the graph is set
         */
        bool HasRoot() const{return has_root_;}
        
        /** 
         * This function returns the size of the grammar count rules by rules
         * @param leaf_bias: how much do leaf nodes count, can be smaller than non-terminal nodes
         * @returns the grammar size of current graph
        */
        double GrammarSize(double leaf_bias);

        /** 
         * This function adds a vertex to the AOG graph
         * @param aog_vertex: the vertex to be added to the AOG graph
         * @throw throw an exception when user tries to add a second root node to the AOG graph
         * @returns an unsigned integer that indicates the id of the newly added vertex(if the vertex alreadly exists, return the existing id).
         */    
        virtual VertexId AddVertex(const std::shared_ptr<AOG_Vertex<StateType, AttributeType> >);

        /**
         * This function deletes a vertex and all related edges from the AOG graph
         * @param vid: the id of the vertex to be deleted
         * @returns if_success true if the vertex is successfully
         */
        virtual bool DeleteVertex(const VertexId);

        /**
         * This function returns all the rules used to construct the AOG graph
         * @returns rules: a vector that contains all the symbolic
         */
        std::vector<Symbolic_Rule<StateType> > GetRules() const;

        /**
         * This function returns all the leaf states in the AOG graph
         * @returns leafStates: a vector that contains all the leaf states
         */
        std::vector<Symbolic_State<StateType> > GetLeafStates() const;

       /**
        * This function returns all the symbolic states in the AOG graph
        * @returns states: a vector that contains all the states in vertices in AOG graph
        */
        std::vector<Symbolic_State<StateType> > GetStates() const;

        /**
         * @param source: the id of the vertex containing the desired state
         * @returns a Symbolic_State object that represents the state inside the given vertex
         */
        Symbolic_State<StateType> GetStateByVertexId(VertexId) const;

        /**
         * Get the root of the AOG
         * @returns vertexId: the vertex id of the root vertex if there exist one, else throw an exception
         */
        VertexId GetRoot() const;

        /**
         * @param state: the symbolic state the user want to query about
         * @returns vid: an unsigned integer that indicates the id of the vertex containing the given symbolic state
         */
        VertexId GetVertexIdByState(const Symbolic_State<StateType> &) const;

        /**
         * This function gets weights of all the outedges from a given source vertex
         * @param source: the id of the source vertex
         * @param is_normalized: a boolean value. If true, return 
         *                     the normalized weights of all the outedges from the target vertex.
         * @returns: weightMap: return an unordered_map that maps all the outedges' target vertices' ids to their corresponding weights
         */
        std::unordered_map<VertexId, double> GetOutEdgeWeights(VertexId, bool) const;

        /**
         * Returns all rules with a certain state as source
         * @param qstate: query state
         * @returns ruleVec: vector of rules with the query state as source, empty if no such rules
        */
        std::vector<Symbolic_Rule<StateType> > GetRulesAsSource(const Symbolic_State<StateType> &) const;

        /**
         * This function get the attribute function for a vertex
         * @param source_id: the id of the query vertex
         * @param neighbors_attrs: the attribute of neighbors' attributes
         * @param args: any other relevant arguments
         */
        std::vector<AttributeType> GetVertexAttributes(VertexId source_id,
            const std::vector<std::vector<AttributeType> >* neighbors_attrs, void *args) const;

        /**
         * Returns all rules with a certain state as one of the targets
         * @param query state
         * @returns vector of rules with the query state as one of the targets, empty if no such rules
        */
        std::vector<Symbolic_Rule<StateType> > GetRulesAsTarget(const Symbolic_State<StateType> &) const;

        /**
         * Returns potential function maps
        */
        const std::unordered_map<std::pair<VertexId, VertexId>,
                                 std::function<double(const std::vector<std::vector<AttributeType> >&, void*)>,
                                 pair_hash>&
        GetBinaryPotentialFuncs() const {return potential_funcs_;};

        /**
         * Returns potential gradient function maps
        */
        const std::unordered_map<std::pair<VertexId, VertexId>,
                                 std::function<void(const std::vector<std::vector<AttributeType> >&, void*, std::vector<double>&)>,
                                 pair_hash>&
        GetBinaryPotentialGradFuncs() const {return potential_grad_funcs_;};

        /**
         * Returns potential function pointer related with the query pair
         * @param query_pair: the query pair of vertices
        */
        const std::unordered_map<std::pair<VertexId, VertexId>,
                             std::function<double(const std::vector<std::vector<AttributeType> >&, void*)>, pair_hash>
        GetBinaryPotentialFunc(const std::pair<VertexId, VertexId>& query_pair) const;

        /**
         * Returns potential gradient function pointer related with the query pair
         * @param query_pair: the query pair of vertices
        */
        const std::unordered_map<std::pair<VertexId, VertexId>,
                             std::function<void(const std::vector<std::vector<AttributeType> >&, void*, std::vector<double>&)>, pair_hash>
        GetBinaryPotentialGradFunc(const std::pair<VertexId, VertexId>& query_pair) const;

        /**
         * Returns potential function parameters related with the query pair
         * @param query_pair: the query pair of vertices
        */
        const std::unordered_map<std::pair<VertexId, VertexId>, void*, pair_hash>
        GetBinaryPotentialParams(const std::pair<VertexId, VertexId>& query_pair) const;

        /**
         * Add a binary potential function to a pair of vertices
         * @param func: the function pointer of the potential function
         * @param vertex_ids: vertex pair related with this potential, order matters
        */
        void SetBinaryPotentialFunc(const std::pair<VertexId, VertexId>& vertex_ids,
                                    const std::function<double(const std::vector<std::vector<AttributeType> >&, void*)> func);

        /**
         * Add a binary gradient function of the potential function to a pair of vertices
         * @param func: the function pointer of the potential function
         * @param vertex_ids: vertex pair related with this potential, order matters
        */
        void SetBinaryPotentialGradFunc(const std::pair<VertexId, VertexId>& vertex_ids,
                                    const std::function<void(const std::vector<std::vector<AttributeType> >&, void*, std::vector<double>&)> func);
        
        /**
         * Set parameters for a binary potential function related with a pair of vertices
         * @param param: the function pointer of the potential function
         * @param vertex_ids: vertex pair related with this potential, order matters
        */
        void SetBinaryPotentialParam(const std::pair<VertexId, VertexId>& vertex_ids, void* param);

        /**
         * Set potential function to a vertex
         * @param func: the function pointer of the potential function
         * @param vertex_id: the vertex related with this potential
        */
        void SetUnaryPotentialFunc(const VertexId vertex_id,
                                   const std::function<double(const std::vector<AttributeType>&, void*)> func);
        
        /**
         * Set gradient function of the potential function to a vertex
         * @param func: the function pointer of the potential function
         * @param vertex_id: the vertex related with this potential
        */
        void SetUnaryPotentialGradFunc(const VertexId vertex_id,
                                       const std::function<void(const std::vector<AttributeType>&, void*, std::vector<double>&)> func);

        /**
         * This function changes a given node to an and/or-node
         * @param source_id: the id of the vertex to be changed
         * @param is_and: true if trying to set the vertex to an and-node
         */
        void SetIsNotOr(VertexId, const bool);

        /**
         * This function set the attribute function for a given node
         * @param source_id: the id of the vertex to be changed
         * @param new_func: the attribute function to be used
         * @return: if set successfully
         */
        bool SetVertexAttributeFunc(VertexId source_id,
                                    const std::function<std::vector<AttributeType>(const AOG_Vertex<StateType, AttributeType>&,
                                                               const std::vector<std::vector<AttributeType> >*, void*)>);

        /**
         * This function set the attribute ranges function for a given node
         * @param source_id: the id of the vertex to be changed
         * @param new_func: the attribute ranges function to be used
         * @param param: the parameter used for this function
         * @return: if set successfully
         */
        bool SetVertexAttributesRangesFunc(VertexId source_id,
            std::function<std::vector<std::vector<AttributeType> >(const AOG_Vertex<StateType, AttributeType>&, void*)> func, void* param);

        /**
         * This function reassigns weights to a given vertex's outedges
         * @param source: the source vertex of all the outedges to be modified
         * @param weight: an unordered_map that contains the mapping from the id of outedges to the weights each outedge need to be set to
         * @returns: true if the weights are set successfully
         */
        bool SetOutEdgeWeights(VertexId, const std::unordered_map<VertexId ,double>&);

        /**
         * This function set the root of the AOG to the vertex holding content root
         * @param root: root state
         */ 
        void SetRoot(Symbolic_State<StateType> root);

        /**
         * This function set the root of the AOG to the vertex holding content root
         * @param self_id: vertex to be set
         * @param ranges: new ranges
         * @return: true if the weights are set successfully
         */
        bool SetVertexAttributesRanges(VertexId source_id, const std::vector<std::vector<AttributeType> >& ranges);

            /**
         * This function normalizes all the weights from a given source vertex
         * @param src_id: the id of the source vertex whose outedges will be normalized
         * @returns an unordered_map that maps ids of all target vertices of outedges to the normalized weights.
         */
        std::unordered_map<VertexId, double> Normalize(VertexId);

        /**
         * This function checks a certain rule for its existence
         * A->BC and A->CB are considered different rules
         * @param rule: the rule to be checked
         * @returns true if the rule already exists in AOG, false otherwise
         */
        bool ExistRule(const Symbolic_Rule<StateType>&);

        /**
         * This function adds a rule to the AOG
         * @param rule: the rule to be added
         * @returns bool: whether the add is successful
         */
        bool AddRule(const Symbolic_Rule<StateType>&);

        /**
         * This function deletes corresponding edges and 
         * vertices in the AOG graph given the rule to be deleted
         * @param rule: the rule to be deleted
         * @returns false if the given rule does not exist, true otherwise
         */
        bool DeleteRule(const Symbolic_Rule<StateType>&);
        
        /**
         * This function merges another AOG into the current AOG
         * by merging all symbolic rules from both side.
         * @param targetAOG: the AOG that wait to be merged.
         * @returns false if two AOGs cannot be merged.
         */
        bool Merge(const AOG<StateType, AttributeType> &);

        /**
         * This Sample() function eliminates all Or-nodes by sampling over edge weights
         * @param root: The vertex to sample from
         * @param sequence: state sequence sampled, pass by reference
         * @param configurations: each element is a configuration (all leaf nodes vertex attributes), pass by reference
         * @returns a vector of vertex ID for the parse graph terminal leaf nodes and
         * sample probability can be obtained with passing by reference
         */
        std::shared_ptr<std::unordered_map<std::pair<VertexId, int>, std::vector<std::pair<VertexId, int> >, pair_hash>> 
        Sample(VertexId root, std::vector<std::pair<VertexId, int> >& sequence,
               std::vector<std::unordered_map<std::pair<VertexId, int>,
                                              std::vector<AttributeType>, pair_hash> >& configurations,
               double& prob, int gibbs_rounds = 100) const;

        /**
         * This GetTopLevelStates() function returns the top most level state(s) of the T-AOG 
         * @returns a vector containing the top most level state(s) of the T-AOG 
         */
        std::vector<Symbolic_State<StateType> > GetTopLevelStates();
        
        /** 
         * This function can be used to remove rules whose source has no parent
         *  The delete is a recursive process
         *  @param state: source state of the rule
         */ 
        void DeleteNoParentRules(const Symbolic_State<StateType>& state);

        /** 
         * Simplify current graph by removing unnecessary dummy nodes introduced
         * in graph editing process without changing the grammar rules.
         * Dummy nodes are internal structure used by this AOG library, doesn't
         * corresponds to any grammar component.
        */
        void TruncateGraph();

        /** 
         * Simplify current graph by compress current grammar 
         * could potentially change the grammar rules, but preserve the semantics
        */
        void SimplifyGraph();

        /**
         * Output a file that can be used for basic AOG visualization
         * @param filename: name of the output file
         * @param dir: name of the output file directory
         * @param truncate: whether remove unnecessary dummy nodes
        */
        void Visualize(std::string dir, std::string filename, bool truncate = true);

        /**
         * Output a file that can be used to reconstruct an AOG
         * @param path: output file path, doesn't include file name, as file name will be learned_tree_appendix.txt
         * @param appendix: the appendix attached in the file name
        */
        void SaveGraph(std::string path, std::string filename);
    };

    template<class StateType, class AttributeType>
    AOG<StateType, AttributeType>::AOG()
    :AOG_Graph<StateType, AttributeType>()
    { this->has_root_ = false; }
    
    //copy constructor for AOG
    template<class StateType, class AttributeType>
    AOG<StateType, AttributeType>::AOG(const AOG &other)
    :AOG_Graph<StateType, AttributeType>(other), all_rules_(other.all_rules_),
    all_leaf_states_(other.all_leaf_states_), state_to_vertex_(other.state_to_vertex_),
    root_id_(other.root_id_), has_root_(other.has_root_)
    {}

    template<class StateType, class AttributeType>
    AOG<StateType, AttributeType>::AOG(const std::vector<Symbolic_Rule<StateType> > &rules)
    :AOG_Graph<StateType, AttributeType>(0, true)
    {
        has_root_ = false;
        for (auto rule: rules)
            this->AddRule(rule);
    }

    template<class StateType, class AttributeType>
    AOG<StateType, AttributeType>::AOG(const std::vector<Symbolic_State<StateType> > &leaf_states)
    :AOG_Graph<StateType, AttributeType>(0, true)
    {
        //construct a root node
        has_root_ = false;
        std::shared_ptr<Symbolic_State<StateType> > root_state =
            std::make_shared<Symbolic_State<StateType> >();
        std::shared_ptr<AOG_Vertex<StateType, AttributeType> > root =
            std::make_shared<AOG_Vertex<StateType, AttributeType> >(*root_state, false, true);
        this->root_id_ = this->AddVertex(root);
        this->state_to_vertex_[*root_state] = this->root_id_;
        has_root_ = true;

        for (unsigned i = 0; i < leaf_states.size(); i++)
        {
            if (!leaf_states[i].GetIsBasic())
            {
                std::cerr << "Non basic states passed in\n";
                throw std::exception();
            }
            //duplicate states in the sequence
            if (this->all_leaf_states_.find(leaf_states[i]) != this->all_leaf_states_.end())
                continue;
            //take in leaf states as content, not an and-node, not a root node
            std::shared_ptr<AOG_Vertex<StateType, AttributeType> > ptr =
                std::make_shared<AOG_Vertex<StateType, AttributeType> >(leaf_states[i], false, false);
            VertexId new_v_id = this->AddVertex(ptr);
            this->all_leaf_states_.insert(leaf_states[i]);
            this->state_to_vertex_[leaf_states[i]] = new_v_id;
        }
    }

    template<class StateType, class AttributeType>
    AOG<StateType, AttributeType>::AOG(const std::string& file_path)
    {
        std::ifstream f(file_path, std::ifstream::in);
        if (f.is_open())
            std::cerr << "reading " << file_path << '\n';
        else
        {
            std::cerr << "file " << file_path <<  " cannot be found" << '\n';
            exit(1);
        }
        
        std::vector<Symbolic_Rule<StateType> > rules;

        std::string s = "";
        int line_num = 0;
        std::unordered_map<Symbolic_Rule<StateType>, double> m;
        while (getline(f, s))
        {
            line_num++;
            std::stringstream ss(s);
            std::string item;
            std::vector<std::string> tokens;
            while (getline(ss, item, ',')) 
                tokens.push_back(item);

            if (tokens.size() < 4)
            {
                std::cerr << "Error: row does not have source id and content\n";
                exit(1);
            }

            int num_res = stoi(tokens[0]);
            double weight = stod(tokens[1]);
            int src_id = stoi(tokens[2]);
            std::string src_ct = tokens[3];
            Symbolic_State<StateType> src(src_id);

            std::vector<Symbolic_State<StateType> > results;
            for (int i = 4; i < tokens.size(); i += 2)
            {
                Symbolic_State<StateType> state;
                int state_id = stoi(tokens[i]);
                if (state_id == -1)
                {
                    std::string state_ct = tokens[i+1];            
                    Symbolic_State<StateType> temp(state_ct, true);
                    state = temp;
                }
                else
                {
                    Symbolic_State<StateType> temp(state_id);
                    state = temp;
                }
                results.push_back(state);
            }
            Symbolic_Rule<StateType> rule(src, results);
            rules.push_back(rule);
            if(weight != 0)
                m[rule] = weight;
        }

        // Find root
        std::unordered_set<Symbolic_State<StateType> > top_level_rules;
        std::vector<Symbolic_State<StateType> > sources;
        std::unordered_set<Symbolic_State<StateType> > results;
        for (Symbolic_Rule<StateType> rule : rules)
        {
            sources.push_back(rule.GetSource());
            results.insert(rule.GetResults().begin(), rule.GetResults().end());
        }
        
        // check which source is not other sources' result
        for (Symbolic_State<StateType> source : sources)
        {
            if (results.find(source) == results.end())
                top_level_rules.insert(source);
        }
        
        if (top_level_rules.size() != 1)
        {
            std::cerr<<"top_level_rule: \n";
            for(const auto & state : top_level_rules)
                std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
            std::cerr << '\n' << "Dangling top level rules beside root\n";
            exit(1);
        }

        std::shared_ptr<AOG_Vertex<StateType, AttributeType> >
            source_ptr(new AOG_Vertex<StateType, AttributeType>(*top_level_rules.begin(), true, true));
        this->AddVertex(source_ptr);

        for (const auto & rule : rules)
            this->AddRule(rule);

        for(const auto& it : m)
        {
            const Symbolic_Rule<StateType>& rule = it.first;
            double weight = it.second;
            VertexId srcId = this->GetVertexIdByState(rule.GetSource());
            std::vector<VertexId> resIds;
            for (const Symbolic_State<StateType>& state : rule.GetResults())
                resIds.push_back(this->GetVertexIdByState(state));
            
            bool found = false;
            VertexId target_dummy;
            for(VertexId dummy : this->ChildrenVertices(srcId))
            {
                std::vector<VertexId> children = this->ChildrenVertices(dummy);
                if(children == resIds)
                {
                    found = true;
                    target_dummy = dummy;
                    break;
                }
            }
            if(!found)
            {
                std::cerr << "rule not found in graph !!" << '\n';
                exit(1);
            }
            auto weights = this->GetOutEdgeWeights(srcId, false);
            weights[target_dummy] = weight;
            this->SetOutEdgeWeights(srcId, weights);
        }
    }

    template<class StateType, class AttributeType>
    std::vector<Symbolic_Rule<StateType> > AOG<StateType, AttributeType>::GetRules() const
    {
        std::vector<Symbolic_Rule<StateType> > all_rules;
        for (auto rule: this->all_rules_)
            all_rules.push_back(rule);
        return all_rules;
    }

    template<class StateType, class AttributeType>
    std::vector<Symbolic_State<StateType> > AOG<StateType, AttributeType>::GetLeafStates() const
    {
        std::vector<Symbolic_State<StateType>> results;
        for (auto it = this->all_leaf_states_.cbegin(); it != this->all_leaf_states_.cend(); it++)
            results.push_back(*it);
        return results;
    }

    template<class StateType, class AttributeType>
    std::vector<Symbolic_State<StateType> > AOG<StateType, AttributeType>::GetStates() const
    {
        std::vector<Symbolic_State<StateType> > all_states(0);
        for (auto state: this->state_to_vertex_)
            all_states.push_back(state.first);
        return all_states;
    }

    template<class StateType, class AttributeType>
    VertexId AOG<StateType, AttributeType>::GetRoot() const 
    {
        if (!this->has_root_){
            std::cout << "The AOG does not have a root\n";
            throw std::exception();
        }
        return this->root_id_;
    }

    template<class StateType, class AttributeType>
    Symbolic_State<StateType> AOG<StateType, AttributeType>::
    GetStateByVertexId(const VertexId vid) const
    {
        if (this->IsValidVertex(vid))
            return this->GetVertexContent(vid)->GetState();
        std::cerr<<"Invalid Vertex Id\n";
        throw std::exception();
    }

    template<class StateType, class AttributeType>
    VertexId AOG<StateType, AttributeType>::
    GetVertexIdByState(const Symbolic_State<StateType> &state) const
    {
        auto iter = this->state_to_vertex_.find(state);
        if (iter == this->state_to_vertex_.end())
        {
            std::cerr << state.GetContent() << " " << state.GetId() << " Entered invalid State\n";
            assert(0);
        }
        return iter->second;
    }

    template<class StateType, class AttributeType>
    std::vector<Symbolic_Rule<StateType> > AOG<StateType, AttributeType>::
    GetRulesAsSource(const Symbolic_State<StateType> & query_state) const
    {
        std::vector<Symbolic_Rule<StateType> > rules;
        for(auto rule: this->all_rules_)
            if(rule.GetSource() == query_state)
                rules.push_back(rule);
        return rules;
    }

    template<class StateType, class AttributeType>
    std::vector<Symbolic_Rule<StateType> > AOG<StateType, AttributeType>::
    GetRulesAsTarget(const Symbolic_State<StateType> & query_state) const
    {
        std::vector<Symbolic_Rule<StateType> > rules;
        for(auto rule: this->all_rules_)
        {
            std::vector<Symbolic_State<StateType> > results = rule.GetResults();
            if(find(results.begin(), results.end(), query_state) != results.end())
                rules.push_back(rule);
        }
        return rules;
    }

    template<class StateType, class AttributeType>
    std::unordered_map<VertexId, double> AOG<StateType, AttributeType>::
    GetOutEdgeWeights(VertexId source, bool is_normalized) const
    {

        //report error if trying to get weights of an and-node's edges
        if(this->GetVertexContent(source)->IsNotOr())
        {
            std::cerr<<"Cannot get weights of edges from an And-node\n";
            throw std::exception();
        }

        std::unordered_map<VertexId, double> weights;
        double total_weights = 0.0;
        for(auto edge:this->OutEdges(source))
        {
            double weight = this->GetEdgeContent(source, edge.second)[0]->GetWeight();
            total_weights += weight;
            weights[edge.second] = weight;
        }

        //if user wants unnormalized weights
        if(!is_normalized)
            return weights;

        double coeff = 1.0/total_weights;
        for(auto &iter : weights)
            iter.second *= coeff;

        return weights;
    }

    template <class StateType, class AttributeType>
    std::vector<AttributeType> AOG<StateType, AttributeType>::
    GetVertexAttributes(VertexId source_id, const std::vector<std::vector<AttributeType> >* neighbors_attrs, void *args) const
    {
        std::vector<AttributeType> results;
        if (!this->IsValidVertex(source_id))
            std::cerr << "Invalid Vertex Id: " << source_id << "\n";
        else if (!(this->GetVertexContent(source_id)->IsNotOr()))
            std::cerr << "Vertex " << source_id << " is or-node.\n";
        else
        {
            if (this->GetVertexContent(source_id)->attributes_ranges_.empty() && (!neighbors_attrs || neighbors_attrs->empty()))
            {
                std::cerr << "[attributes are either sampled from self attribute's range or calculated with neighbor's attributes,"
                             "if sampled, initialize self attributes' range and pass in null neighbors' attribute pointer,"
                             "otherwise pass in valid neighbors' attributes]\n";
                assert(0);
            }
            return this->GetVertexContent(source_id)->GetAttributes(neighbors_attrs, args);
        }
        return results;
    }
    
    template <class StateType, class AttributeType>
    const std::unordered_map<std::pair<VertexId, VertexId>,
                             std::function<double(const std::vector<std::vector<AttributeType> >&, void*)>, pair_hash>
    AOG<StateType, AttributeType>::
    GetBinaryPotentialFunc(const std::pair<VertexId, VertexId>& query_pair) const
    {
        std::unordered_map<std::pair<VertexId, VertexId>,
                           std::function<double(const std::vector<std::vector<AttributeType> >&, void*)>, pair_hash> results;
        VertexId vid1 = query_pair.first;
        VertexId vid2 = query_pair.second;
        if (!(this->IsValidVertex(vid1) && this->IsValidVertex(vid2)))
            std::cerr << "Vertex:" << vid1 << " or vertex: " << vid2 << " is invalid\n";
        else
        {
            std::vector<VertexId> element1;
            std::vector<VertexId> element2;
            if(this->GetVertexContent(vid1)->IsNotOr())
                element1.push_back(vid1);
            else
                for (auto vid : this->ChildrenVertices(vid1))
                    element1.push_back(vid);
            if(this->GetVertexContent(vid2)->IsNotOr())
                element2.push_back(vid2);
            else
                for (auto vid : this->ChildrenVertices(vid2))
                    element2.push_back(vid);
            for (auto e1: element1)
                for(auto e2: element2)
                    if (this->potential_funcs_.count(std::make_pair(e1, e2)))
                        results[std::make_pair(e1, e2)] = (this->potential_funcs_.at(std::make_pair(e1, e2)));
        }
        return results;
    }

    template <class StateType, class AttributeType>
    const std::unordered_map<std::pair<VertexId, VertexId>,
                             std::function<void(const std::vector<std::vector<AttributeType> >&, void*, std::vector<double>&)>, pair_hash>
    AOG<StateType, AttributeType>::
    GetBinaryPotentialGradFunc(const std::pair<VertexId, VertexId>& query_pair) const
    {
        std::unordered_map<std::pair<VertexId, VertexId>,
                           std::function<void(const std::vector<std::vector<AttributeType> >&, void*, std::vector<double>&)>, pair_hash> results;
        VertexId vid1 = query_pair.first;
        VertexId vid2 = query_pair.second;
        if (!(this->IsValidVertex(vid1) && this->IsValidVertex(vid2)))
            std::cerr << "Vertex:" << vid1 << " or vertex: " << vid2 << " is invalid\n";
        else
        {
            std::vector<VertexId> element1;
            std::vector<VertexId> element2;
            if(this->GetVertexContent(vid1)->IsNotOr())
                element1.push_back(vid1);
            else
                for (auto vid : this->ChildrenVertices(vid1))
                    element1.push_back(vid);
            if(this->GetVertexContent(vid2)->IsNotOr())
                element2.push_back(vid2);
            else
                for (auto vid : this->ChildrenVertices(vid2))
                    element2.push_back(vid);
            for (auto e1: element1)
                for(auto e2: element2)
                    if (this->potential_funcs_.count(std::make_pair(e1, e2)))
                        results[std::make_pair(e1, e2)] = (this->potential_grad_funcs_.at(std::make_pair(e1, e2)));
        }
        return results;
    }

    template <class StateType, class AttributeType>
    const std::unordered_map<std::pair<VertexId, VertexId>, void*, pair_hash> AOG<StateType, AttributeType>::
    GetBinaryPotentialParams(const std::pair<VertexId, VertexId>& query_pair) const
    {
        std::unordered_map<std::pair<VertexId, VertexId>, void*, pair_hash> results;
        VertexId vid1 = query_pair.first;
        VertexId vid2 = query_pair.second;
        if (!(this->IsValidVertex(vid1) && this->IsValidVertex(vid2)))
            std::cerr << "Vertex:" << vid1 << " or vertex: " << vid2 << " is invalid\n";
        else
        {
            std::vector<VertexId> element1;
            std::vector<VertexId> element2;
            if(this->GetVertexContent(vid1)->IsNotOr())
                element1.push_back(vid1);
            else
                for (auto vid : this->ChildrenVertices(vid1))
                    element1.push_back(vid);
            if(this->GetVertexContent(vid2)->IsNotOr())
                element2.push_back(vid2);
            else
                for (auto vid : this->ChildrenVertices(vid2))
                    element2.push_back(vid);
            for (auto e1: element1)
                for(auto e2: element2)
                    if (this->potential_params_.count(std::make_pair(e1, e2)))
                        results[std::make_pair(e1, e2)] = (this->potential_params_.at(std::make_pair(e1, e2)));
        }
        return results;
    }

    template <class StateType, class AttributeType>
    void AOG<StateType, AttributeType>::
    SetBinaryPotentialFunc(const std::pair<VertexId, VertexId>& query_pair,
                           const std::function<double(const std::vector<std::vector<AttributeType> >&, void*)> func)
    {
        VertexId vid1 = query_pair.first;
        VertexId vid2 = query_pair.second;
        if (!(this->IsValidVertex(vid1) && this->IsValidVertex(vid2)))
        {
            std::cerr << "Vertex:" << vid1 << " or vertex: " << vid2 << " is invalid\n";
            assert(0);
        }
        else
        {
            std::vector<VertexId> element1;
            std::vector<VertexId> element2;
            if(this->GetVertexContent(vid1)->IsNotOr())
                element1.push_back(vid1);
            else
                for (auto cid : this->ChildrenVertices(vid1))
                    element1.push_back(cid);
            if(this->GetVertexContent(vid2)->IsNotOr())
                element2.push_back(vid2);
            else
                for (auto cid : this->ChildrenVertices(vid2))
                    element2.push_back(cid);
            for (auto e1: element1)
                for(auto e2: element2)
                {
                    this->potential_funcs_[std::make_pair(e1, e2)] = func;
                    this->potential_params_[std::make_pair(e1, e2)] = 0;
                }
        }
    }

    template <class StateType, class AttributeType>
    void AOG<StateType, AttributeType>::
    SetBinaryPotentialGradFunc(const std::pair<VertexId, VertexId>& query_pair,
                               const std::function<void(const std::vector<std::vector<AttributeType> >&, void*, std::vector<double>&)> func)
    {
        VertexId vid1 = query_pair.first;
        VertexId vid2 = query_pair.second;
        if (!(this->IsValidVertex(vid1) && this->IsValidVertex(vid2)))
        {
            std::cerr << "Vertex:" << vid1 << " or vertex: " << vid2 << " is invalid\n";
            assert(0);
        }
        else
        {
            std::vector<VertexId> element1;
            std::vector<VertexId> element2;
            if(this->GetVertexContent(vid1)->IsNotOr())
                element1.push_back(vid1);
            else
                for (auto cid : this->ChildrenVertices(vid1))
                    element1.push_back(cid);
            if(this->GetVertexContent(vid2)->IsNotOr())
                element2.push_back(vid2);
            else
                for (auto cid : this->ChildrenVertices(vid2))
                    element2.push_back(cid);
            for (auto e1: element1)
                for(auto e2: element2)
                {
                    this->potential_grad_funcs_[std::make_pair(e1, e2)] = func;
                }
        }
    }

    template <class StateType, class AttributeType>
    void AOG<StateType, AttributeType>::
    SetBinaryPotentialParam(const std::pair<VertexId, VertexId>& query_pair, void* params)
    {
        VertexId vid1 = query_pair.first;
        VertexId vid2 = query_pair.second;
        if (!(this->IsValidVertex(vid1) && this->IsValidVertex(vid2)))
        {
            std::cerr << "Vertex:" << vid1 << " or vertex: " << vid2 << " is invalid\n";
            assert(0);
        }
        else
        {
            std::vector<VertexId> element1;
            std::vector<VertexId> element2;
            if(this->GetVertexContent(vid1)->IsNotOr())
                element1.push_back(vid1);
            else
                for (auto cid : this->ChildrenVertices(vid1))
                    element1.push_back(cid);
            if(this->GetVertexContent(vid2)->IsNotOr())
                element2.push_back(vid2);
            else
                for (auto cid : this->ChildrenVertices(vid2))
                    element2.push_back(cid);
            for (auto e1: element1)
                for(auto e2: element2)
                {
                    this->potential_params_[std::make_pair(e1, e2)] = params;
                }
        }
    }

    template <class StateType, class AttributeType>
    void AOG<StateType, AttributeType>::
    SetUnaryPotentialFunc(const VertexId vertex_id, const std::function<double(const std::vector<AttributeType>&, void*)> func)
    {
        if (!(this->IsValidVertex(vertex_id)))
        {
            std::cerr << "Vertex:" << vertex_id << " is invalid\n";
            assert(0);
        }
        else
        {
            if(this->GetVertexContent(vertex_id)->IsNotOr())
                this->GetVertexContent(vertex_id)->SetPotentialFunc(func);
            else
                for (auto cid : this->ChildrenVertices(vertex_id))
                    this->GetVertexContent(cid)->SetPotentialFunc(func);
        }
    }

    template <class StateType, class AttributeType>
    void AOG<StateType, AttributeType>::
    SetUnaryPotentialGradFunc(const VertexId vertex_id,
                              const std::function<void(const std::vector<AttributeType>&, void*, std::vector<double>&)> func)
    {
        if (!(this->IsValidVertex(vertex_id)))
        {
            std::cerr << "Vertex:" << vertex_id << " is invalid\n";
            assert(0);
        }
        else
        {
            if(this->GetVertexContent(vertex_id)->IsNotOr())
                this->GetVertexContent(vertex_id)->SetPotentialGradFunc(func);
            else
                for (auto cid : this->ChildrenVertices(vertex_id))
                    this->GetVertexContent(cid)->SetPotentialGradFunc(func);
        }
    }

    //T-AOG version add vertex, avoid adding duplicate vertices if they have same state
    template <class StateType, class AttributeType>
    VertexId AOG<StateType, AttributeType>::
        AddVertex(const std::shared_ptr<AOG_Vertex<StateType, AttributeType>> aog_vertex)
    {
        Symbolic_State<StateType> state = aog_vertex->GetState();
        bool is_and = aog_vertex->IsNotOr();
        bool is_root = aog_vertex->IsRoot();
        //check if user wants to add a second root
        if(is_root && has_root_)
        {
            std::cerr << "Cannot add another root vertex.\n";
            throw std::exception();
        }
        //if the symbolic state in the given vertex does not exist, add a new vertex
        auto iter = this->state_to_vertex_.find(state);
        
        if (iter == this->state_to_vertex_.end())
        {
            std::shared_ptr<AOG_Vertex<StateType, AttributeType> >
                ptr(new AOG_Vertex<StateType, AttributeType>(state, is_and, is_root));
            ptr->attributes_ranges_ = aog_vertex->attributes_ranges_;
            ptr->attribute_func_ = aog_vertex->attribute_func_;
            ptr->potential_func_ = aog_vertex->potential_func_;
            ptr->potential_params_ = aog_vertex->potential_params_;
            VertexId new_id = AOG_Graph<StateType, AttributeType>::AddVertex(ptr);
            this->state_to_vertex_[state] = new_id;
            if (state.GetIsBasic())
            {
                assert(state.GetContent().length() != 0);
                this->all_leaf_states_.insert(state);
            }
            
            //if the vertex added want to be a root, add this vertex as root
            if(is_root && !(this->has_root_))
            {
                this->has_root_ = true;
                this->root_id_ = new_id;
            }
            return new_id;
        }

        //else return the existing vertex's id
        if(is_root && !(this->has_root_))
        {
            this->has_root_ = true;
            this->root_id_ = iter->second;
        }
        return iter->second;
    }

    template<class StateType, class AttributeType>
    bool AOG<StateType, AttributeType>::DeleteVertex(const VertexId vid)
    {
        //update the data structures that keep track of all states and leaf states
        if (this->IsValidVertex(vid))
        {
            std::shared_ptr<AOG_Vertex<StateType, AttributeType> > aog_vertex = this->GetVertexContent(vid);
            Symbolic_State<StateType> state = aog_vertex->GetState();
            if (state.GetIsBasic())
                this->all_leaf_states_.erase(state);
            if(this->state_to_vertex_.at(state) == vid)
                this->state_to_vertex_.erase(state);
        }

        //Invalid situation handled in base class
        return AOG_Graph<StateType, AttributeType>::DeleteVertex(vid);
    }


    template<class StateType, class AttributeType>
    bool AOG<StateType, AttributeType>::AddEdge(const VertexId source, const VertexId target,
                                                const std::shared_ptr<AOG_Edge> aog_edge,
                                                bool multi_edge, int pos)
    {
        double weight = aog_edge->GetWeight();
        //bound checking
        if(weight < 0)
        {
            std::cerr << "Weight cannot be less than 0!\n";
            return false;
        }
        if(this->IsValidVertex(source))
        {
            //source of the edge to add should not be a leaf state
            Symbolic_State<StateType> src_state = GetStateByVertexId(source);
            if(src_state.GetIsBasic())
            {
                std::cerr<<"cannot add edge starting from leaf state\n";
                return false;
            }
        }

        //invalid situations handled in base class
        std::shared_ptr<AOG_Edge> ptr = std::make_shared<AOG_Edge>(weight);
        return AOG_Graph<StateType, AttributeType>::AddEdge(source, target, ptr, multi_edge, pos);
    }

    template <class StateType, class AttributeType>
    void AOG<StateType, AttributeType>::SetRoot(Symbolic_State<StateType> root)
    {
        this->has_root_ = true;
        this->root_id_ = this->GetVertexIdByState(root);
        this->GetVertexContent(this->root_id_)->SetIsRoot(true);
    }

    template<class StateType, class AttributeType>
    bool AOG<StateType, AttributeType>::
    SetOutEdgeWeights(VertexId source, const std::unordered_map<VertexId, double> & weights)
    {

        if(this->GetVertexContent(source)->IsNotOr())
        {
            std::cerr<<"Cannot get weights of edges from an And-node\n";
            return false;
        }
        //check if the number of weights inputs is the same as the number of out edges
        if(weights.size() != this->OutEdges(source).size())
        {
            std::cerr<<"The number of weights inputs is inconsistent with the number of out edges.\n";
            return false;
        }
        //update weights
        for(auto iter:weights)
            this->GetEdgeContent(source, iter.first)[0]->SetWeight(iter.second);

        return true;
    }

    template<class StateType, class AttributeType>
    void AOG<StateType, AttributeType>::SetIsNotOr(VertexId source_id, const bool is_and)
    {
        if (this->IsValidVertex(source_id))
            this->GetVertexContent(source_id)->SetIsNotOr(is_and);
        else
            std::cerr << "Vertex" << source_id << " doesn't exist.\n";

        return;
    }

    template<class StateType, class AttributeType>
    bool AOG<StateType, AttributeType>::
    SetVertexAttributeFunc(VertexId source_id,
                           const std::function<std::vector<AttributeType>(const AOG_Vertex<StateType, AttributeType>&,
                                                      const std::vector<std::vector<AttributeType> >*, void*)> func)
    {
        if(this->IsValidVertex(source_id) && this->GetVertexContent(source_id)->IsNotOr())
            this->GetVertexContent(source_id)->SetAttributeFunc(func);
        else if(!this->IsValidVertex(source_id))
        {
            std::cerr << "Vertex" << source_id << " doesn't exist.\n";
            return false;
        }
        else
        {
            for (auto cid : this->ChildrenVertices(source_id))
                this->GetVertexContent(cid)->SetAttributeFunc(func);
        }
        return true;
    }


    template <class StateType, class AttributeType>
    bool AOG<StateType, AttributeType>::SetVertexAttributesRangesFunc(VertexId source_id,
        std::function<std::vector<std::vector<AttributeType> >(const AOG_Vertex<StateType, AttributeType>&, void*)> func, void* param)
    {
        if (this->IsValidVertex(source_id) && this->GetVertexContent(source_id)->IsNotOr())
            this->GetVertexContent(source_id)->SetAttributesRangesFunc(func, param);
        else if (!this->IsValidVertex(source_id))
        {
            std::cerr << "Vertex" << source_id << " doesn't exist.\n";
            return false;
        }
        else
        {
            for (auto cid : this->ChildrenVertices(source_id))
                this->GetVertexContent(cid)->SetAttributesRangesFunc(func, param);
        }
        return true;
    }

    template <class StateType, class AttributeType>
    bool AOG<StateType, AttributeType>::SetVertexAttributesRanges(VertexId source_id,
        const std::vector<std::vector<AttributeType> >& ranges)
    {
        if (this->IsValidVertex(source_id) && this->GetVertexContent(source_id)->IsNotOr())
            this->GetVertexContent(source_id)->SetAttributesRanges(ranges);
        else if (!this->IsValidVertex(source_id))
        {
            std::cerr << "Vertex" << source_id << " doesn't exist.\n";
            return false;
        }
        else
        {
            for (auto cid : this->ChildrenVertices(source_id))
                this->GetVertexContent(cid)->SetAttributesRanges(ranges);
        }
        return true;
    }

    template<class StateType, class AttributeType>
    std::unordered_map<VertexId,double> AOG<StateType, AttributeType>::Normalize(VertexId src_id)
    {
        //if it is an And-node, do nothing
        if(this->GetVertexContent(src_id)->IsNotOr())
        {
            std::cerr << "Cannot normalize edges from And-node\n";
            throw std::exception();
        }

        //get normalized weights
        std::unordered_map<VertexId,double> weights = this->GetOutEdgeWeights(src_id,true);
        //set weights
        this->SetOutEdgeWeights(src_id,weights);
        return weights;
    }
    
    template<class StateType, class AttributeType>
    bool AOG<StateType, AttributeType>::ExistRule(const Symbolic_Rule<StateType>& rule)
    {
        return this->all_rules_.find(rule) != this->all_rules_.end();
    }

    template<class StateType, class AttributeType>
    bool AOG<StateType, AttributeType>::AddRule(const Symbolic_Rule<StateType>& rule)
    {
        // construct an edge with weight 1.0
        std::shared_ptr<AOG_Edge> edge = std::make_shared<AOG_Edge>();

        // check if the rule already exists
        if (this->ExistRule(rule))        
            return false;

        // get the source and result states from each rule and construct vertices for each of them.
        std::shared_ptr<Symbolic_State<StateType> > source_state = 
            std::make_shared<Symbolic_State<StateType> >(rule.GetSource());
        std::shared_ptr<AOG_Vertex<StateType, AttributeType> > source = 
            std::make_shared<AOG_Vertex<StateType, AttributeType> >(*source_state, true, false);
        
        // std::cerr<<"this rule is added from state id: "<<(*source_state).GetId()<<std::endl;
        VertexId src_id = this->AddVertex(source);               
        this->state_to_vertex_[*source_state] = src_id;
    
        std::vector<std::shared_ptr<Symbolic_State<StateType> > > result_states;
        std::vector<VertexId> rs_vtxs_id;
        auto results = rule.GetResults();

        for(auto result : results)
        {
            // std::cerr<<"this state is added to state id: "<<result.GetId();
            result_states.push_back(std::make_shared<Symbolic_State<StateType> >(result));
            std::shared_ptr<AOG_Vertex<StateType, AttributeType> > rs_vtx = 
                std::make_shared<AOG_Vertex<StateType, AttributeType> >(result, true, false);
            VertexId rs_id = this->AddVertex(rs_vtx);
            rs_vtxs_id.push_back(rs_id);
            if(result.GetIsBasic())
            {
                assert(result.GetContent().length() != 0);
                this->all_leaf_states_.insert(result);
            }
            this->state_to_vertex_[result] = rs_id;
        }
                       
        std::vector<VertexId> children = AOG_Graph<StateType, AttributeType>::ChildrenVertices(src_id);


        //create a dummy And-node
        std::shared_ptr<Symbolic_State<StateType> > dummy_state = 
            std::make_shared<Symbolic_State<StateType> >(*source_state);
        std::shared_ptr<AOG_Vertex<StateType, AttributeType> > dummy = 
            std::make_shared<AOG_Vertex<StateType, AttributeType> >(*dummy_state, true, false);
        dummy->attributes_ranges_ = this->GetVertexContent(src_id)->attributes_ranges_;
        dummy->attribute_func_ = this->GetVertexContent(src_id)->attribute_func_;
        dummy->potential_func_ = this->GetVertexContent(src_id)->potential_func_;
        dummy->potential_params_ = this->GetVertexContent(src_id)->potential_params_;
        // if src is already a parent
        if (children.size())
        {
            // Add the dummy And-node for the new rule
            VertexId dummy_id = AOG_Graph<StateType, AttributeType>::AddVertex(dummy);

            // if src is an And-node, it has no dummy nodes yet
            if (this->GetVertexContent(src_id)->IsNotOr())
            {
                // create an Or-node new parent between src and its parent
                std::shared_ptr<Symbolic_State<StateType> > new_par_state =
                    std::make_shared<Symbolic_State<StateType> >(*source_state);
                std::shared_ptr<AOG_Vertex<StateType, AttributeType> > new_par =
                    std::make_shared<AOG_Vertex<StateType, AttributeType> >(*new_par_state, false, false);
                // no need to copy attribute related member variable as new parent is an or-noce
                //delete the old state and add the new state
                this->state_to_vertex_.erase(*source_state);
                VertexId new_par_id = this->AddVertex(new_par);
                this->state_to_vertex_[*new_par_state] = new_par_id;
                //change the root to the newly created node
                if(this->has_root_ && src_id == this->GetRoot())
                {
                    this->root_id_ = new_par_id;
                    this->GetVertexContent(src_id)->SetIsRoot(false);
                    this->GetVertexContent(new_par_id)->SetIsRoot(true);

                }

                // if src has parents, connect parents to new_par
                // delete parents to src
                for (auto parent_id : this->ParentsVertices(src_id)) 
                {                        
                    //find the source node's positions in all rules whose targets contains source
                    std::vector<VertexId> all_children = this->ChildrenVertices(parent_id);                
                    std::vector<int> result_pos;
                    auto result = std::find(all_children.begin(),all_children.end(),src_id);
                    while(result != all_children.end())
                    {
                        result_pos.push_back(result - all_children.begin());
                        result = std::find(result + 1,all_children.end(),src_id);
                    }
                    
                    //delete all edges that have given source and target
                    this->DeleteEdge(parent_id, src_id);

                    //add all edges from parent to newly added state
                    for(int pos : result_pos)
                        this->AddEdge(parent_id, new_par_id, edge, true, pos);

                }
                // connect new_par to src and dummy
                this->AddEdge(new_par_id, src_id, edge);
                this->AddEdge(new_par_id, dummy_id, edge);
                // connect dummy to all results
                for(auto rs_id : rs_vtxs_id)
                    this->AddEdge(dummy_id, rs_id, edge);

            }
            // if src is an Or-node, simply add another dummy node
            else
            {
                // connect src to dummy
                this->AddEdge(src_id, dummy_id, edge);
                // connect dummy to all results
                for(auto rs_id : rs_vtxs_id)
                    this->AddEdge(dummy_id, rs_id, edge);
            }
        }

        // if src is not yet a parent, and src is root,
        // add a dummy node below, set root to or_node, then add children under dummy  
        else if(!children.size() && this->GetVertexContent(src_id)->IsRoot()) 
        {
            VertexId dummy_id = AOG_Graph<StateType, AttributeType>::AddVertex(dummy);
            this->AddEdge(src_id, dummy_id, edge);
           
            for(auto rs_id : rs_vtxs_id)
                this->AddEdge(dummy_id, rs_id, edge);
            
            //change root to or-node
            this->GetVertexContent(src_id)->SetIsNotOr(false);

        }

        // if src is not yet a parent, and src is not root, just add edges to its children
        else
        {
            // std::cerr<<"the edge is added from: "<<src_id<<std::endl;
            for(auto rs_id : rs_vtxs_id)
            {
                // std::cerr<<"the edge is added to : "<<rs_id<<std::endl;
                this->AddEdge(src_id, rs_id, edge);
            }
            
        }

        this->all_rules_.insert(rule);
        return true;
    }
    
    template<class StateType, class AttributeType>
    bool AOG<StateType, AttributeType>::DeleteRule(const Symbolic_Rule<StateType>& rule)
    {
      //delete the rule in all_rules_
        if(!all_rules_.erase(rule))
        {
            //std::cerr<<"The rule to be deleted does not exist in the AOG!\n";
            return false;
        }
    
        //get vertices ids
        VertexId src_vtx = GetVertexIdByState(rule.GetSource());
        auto results = rule.GetResults();
        std:: vector<VertexId> end_vtxs;
        for (auto result : results)
            end_vtxs.push_back(GetVertexIdByState(result));
    
        size_t end_vtxs_size = end_vtxs.size();
    

        //if the source node is an and-node
        if(this->GetVertexContent(src_vtx)->IsNotOr())
        {
            //delete the edges constructed from the rule
            for(auto end_vtx : end_vtxs)
                this->DeleteEdge(src_vtx,end_vtx);
        }

        //if source state is an or-node
        else
        {
            //find the dummy and delete it
            auto dummys = this->ChildrenVertices(src_vtx);
            bool found_rule;
            for(auto dummy : dummys)
            {
                found_rule = true;
                auto dummy_children = this->ChildrenVertices(dummy);
        
                if(dummy_children.size() != end_vtxs_size)
                    continue;
                //if the child is the dummy node that associated with the rule, delete the dummy
                for(size_t i = 0; i< end_vtxs_size;i++)
                    if(dummy_children[i] != end_vtxs[i])
                    {
                        found_rule = false;
                        break;
                    }
                if(found_rule)
                {
                    AOG_Graph<StateType, AttributeType>::DeleteVertex(dummy);
                    break;
                }
            }

            //if there is only one dummy left, and the source is not root, change the source vertex to an and-node
            dummys = this->ChildrenVertices(src_vtx);
            if(dummys.size() == 1 && (!this->GetVertexContent(src_vtx)->IsRoot()))
            {
                std::vector<VertexId> dummy_children = this -> ChildrenVertices(dummys[0]);
                for(auto dummy_child : dummy_children)
                    AddEdge(src_vtx,dummy_child,std::make_shared<AOG_Edge>());
                //delete the dummy node
                AOG_Graph<StateType, AttributeType>::DeleteVertex(dummys[0]);
        
                //change the property of the source vertex to and-node
                this->GetVertexContent(src_vtx)->SetIsNotOr(true);
            }
        }
        return true;
    }

    template <class StateType, class AttributeType>
    std::shared_ptr<std::unordered_map<std::pair<VertexId, int>, std::vector<std::pair<VertexId, int> >, pair_hash> >
    AOG<StateType, AttributeType>::Sample(VertexId root, std::vector<std::pair<VertexId, int> >& res,
                                          std::vector<std::unordered_map<std::pair<VertexId, int>,
                                                                         std::vector<AttributeType>, pair_hash> >& configurations,
                                          double &prob, int gibbs_rounds) const
    {
        prob = 1;
        //if the node to sample from is a dummy node, start sampling from its parent
        auto root_parent = this->ParentsVertices(root);
        if(!root_parent.empty() && !this->GetVertexContent(root_parent[0])->IsNotOr())
            root = root_parent[0];
        // the resultant parse graph that contains the depth-nodes pairs
        std::shared_ptr<std::unordered_map<std::pair<VertexId, int>, std::vector<std::pair<VertexId, int> >, pair_hash> > parse_tree =
            std::make_shared<std::unordered_map<std::pair<VertexId, int>, std::vector<std::pair<VertexId, int> >, pair_hash> >();
        std::unordered_map<VertexId, int> occurance_count;
        
        std::stack<std::pair<VertexId, int> > stack;
        std::deque<std::pair<VertexId, int> > partial_order; //partial order of the non-terminal-vertices
        // one binary potential function can be triggered by multiple pairs
        std::unordered_map<std::pair<VertexId, VertexId>, std::vector<std::pair<std::pair<VertexId, int>, std::pair<VertexId, int> > >, pair_hash>
            included_binary_potential_funcs;
        std::vector<std::pair<VertexId, int> > included_unary_potential_funcs;
        // will be used to obtain a seed for the random number engine
        std::random_device rd;
        // standard mersenne_twister_engine seeded with rd()
        std::mt19937 gen(rd());

        stack.push(std::make_pair(root,0));
        occurance_count.emplace(root, 0);
        occurance_count[root]++;
        auto update_included_binary_potential = [&](const std::pair<VertexId, int>& new_vertex_pair)
        {
            VertexId new_vertex_id = new_vertex_pair.first;
            for (auto it = occurance_count.cbegin(); it != occurance_count.cend(); it++)
            {
                std::pair<VertexId, VertexId> pair1 = std::make_pair(new_vertex_id, it->first);
                std::pair<VertexId, VertexId> pair2 = std::make_pair(it->first, new_vertex_id);
                if (this->potential_funcs_.count(pair1))
                {
                    for (int oc = 0; oc < it->second; oc++)
                    {
                        included_binary_potential_funcs.emplace(pair1,
                            std::vector<std::pair<std::pair<VertexId, int>, std::pair<VertexId, int> > >());
                        included_binary_potential_funcs.at(pair1).push_back(
                            std::make_pair(new_vertex_pair, std::make_pair(it->first, oc)));
                    }
                }
                else if (this->potential_funcs_.count(pair2))
                {
                    for (int oc = 0; oc < it->second; oc++)
                    {
                        included_binary_potential_funcs.emplace(pair2,
                            std::vector<std::pair<std::pair<VertexId, int>, std::pair<VertexId, int> > >());
                        included_binary_potential_funcs.at(pair2).push_back(
                            std::make_pair(std::make_pair(it->first, oc), new_vertex_pair));
                    }
                }
            }
        };
        while (!stack.empty())
        {
            std::pair<VertexId, int> parent_pair = stack.top();
            VertexId parent_id = parent_pair.first;
            int parent_id_count = parent_pair.second;
            stack.pop();
            std::vector<VertexId> children = this->ChildrenVertices(parent_id);
            if (children.size() == 0 && this->GetVertexContent(parent_pair.first)->GetState().GetIsBasic())
            {
                res.push_back(parent_pair);
            }
            // if the node is an And-node, keep all its children in sample_graph
            else if (AOG_Graph<StateType, AttributeType>::GetVertexContent(parent_id)->IsNotOr())
            {
                partial_order.push_front(parent_pair);
                for (int i = children.size() - 1; i >= 0; i--)
                {
                    occurance_count.emplace(children[i], 0);
                    std::pair<VertexId, int> child_pair = std::make_pair(children[i], occurance_count.at(children[i])++);
                    stack.push(child_pair);
                    parse_tree->emplace(parent_pair, std::vector<std::pair<VertexId, int> >());
                    parse_tree->at(parent_pair).insert(parse_tree->at(parent_pair).begin(), child_pair);
                    if (this->GetVertexContent(children[i])->IsNotOr())
                    {
                        // add unary potential function
                        if (this->GetVertexContent(children[i])->potential_func_)
                            included_unary_potential_funcs.push_back(child_pair);
                        //add binary potential function
                        update_included_binary_potential(child_pair);
                    }
                }
            }
            // if the node is an Or-node
            else
            {
                std::unordered_map<VertexId, double> edge_weights = this->GetOutEdgeWeights(parent_id, true);
                
                std::vector<VertexId> vertex;
                std::vector<double> weight;

                for (auto w:edge_weights)
                {
                    vertex.push_back(w.first);
                    weight.push_back(w.second);
                }

                // choose a child under the weight distribution and keep it in sample_graph
                std::discrete_distribution<> dis(weight.begin(), weight.end());
                unsigned index = dis(gen);
                VertexId sample = vertex[index];
                prob *= weight[index];
                
                occurance_count.emplace(sample, 0);
                std::pair<VertexId, int> sample_pair = std::make_pair(sample, occurance_count.at(sample)++);
                stack.push(sample_pair);
                
                //substitute previous or-node with its child
                for(auto it = parse_tree->begin(); it != parse_tree->end(); it++)
                    for(auto& child_nodes: it->second)
                        if(child_nodes == parent_pair)
                            child_nodes = sample_pair;
                if (this->GetVertexContent(sample)->potential_func_)
                            included_unary_potential_funcs.push_back(sample_pair);
                update_included_binary_potential(sample_pair);
            }
        }

        std::unordered_map<std::pair<VertexId, int>, std::vector<AttributeType>, pair_hash> configuration;
        int initial_terminal_arg[1] = {-1};
        for (auto re: res) 
            if (this->GetVertexContent(re.first)->attribute_func_)
                configuration[re] = this->GetVertexAttributes(re.first, 0, initial_terminal_arg);
        
        auto bottom_up = [&](std::unordered_map<std::pair<VertexId, int>, std::vector<AttributeType>, pair_hash>& configure)
        {
            for (auto po: partial_order)
            {
                if (this->GetVertexContent(po.first)->attribute_func_)
                {
                    std::vector<std::vector<AttributeType> > nei_attr;
                    for (auto nei: parse_tree->at(po))
                    {
                        if (this->GetVertexContent(nei.first)->attribute_func_)
                            nei_attr.push_back(configure.at(nei));
                    }
                    configure[po] = this->GetVertexAttributes(po.first, &nei_attr, 0);
                }
            }
        };

        auto calculate_score = [&](std::unordered_map<std::pair<VertexId, int>, std::vector<AttributeType>, pair_hash>& configure)
        {
            double score = 0;
            for (auto iup: included_unary_potential_funcs)
                score += this->GetVertexContent(iup.first)->
                    GetPotentialFunc()(configure[iup], this->GetVertexContent(iup.first)->GetPotentialParam());
            for (auto ibp_it = included_binary_potential_funcs.cbegin(); ibp_it != included_binary_potential_funcs.cend(); ibp_it++)
                for (auto couple: ibp_it->second)
                {
                    if (configure.count(couple.first) == 0 || configure.count(couple.second) == 0)
                    {
                        std::cerr << "cannot set potential function over vertices without attributes\n";
                        assert(0);
                    }
                    score += this->GetBinaryPotentialFuncs().at(ibp_it->first)
                        ({configure.at(couple.first), configure.at(couple.second)},
                         this->potential_params_.at(std::make_pair(couple.first.first, couple.second.first)));
                }
            return score;
        };

        bottom_up(configuration);
        configurations.push_back(configuration);

        auto neg_softmax = [](const std::vector<double>& logits)
        {
            std::vector<double> probs(logits.size(), 0);
            double total_prob = 0;
            double mean_l = 0;
            for(auto logit: logits)
                mean_l += logit;
            mean_l /= logits.size();
            for (int i = 0; i < logits.size(); i++)
            {
                probs[i] = exp(-1 * (logits[i] - mean_l));
                total_prob += probs[i];
            }
            for (int i = 0; i < logits.size(); i++)
                probs[i] /= total_prob;
            return probs;
        };
        auto start_time = Clock::now();
        for (int gr = 0; gr < gibbs_rounds; gr++)
        {
            std::unordered_map<std::pair<VertexId, int>, std::vector<AttributeType>, pair_hash>
            current_configuration = configurations.back();
            // gibbs only change attributes of the leaf nodes, non-terminal nodes are calculated
            for (auto re: res)
            {
                std::vector<std::vector<AttributeType> > ranges = this->GetVertexContent(re.first)->attributes_ranges_;
                //for each dimension of the attribute
                for (int a_pos = 0; a_pos < ranges.size(); a_pos++)
                {
                    int num_possibles = ranges[a_pos].size();
                    //for each element in that attribute dimension range
                    int gibbs_args[2] = {a_pos, -1};
                    std::vector<std::vector<AttributeType> > current_attribute(1, {current_configuration.at(re)});
                    std::vector<double> energies;
                    for (int a_select = 0; a_select < num_possibles; a_select++)
                    {
                        gibbs_args[1] = a_select;
                        current_configuration.at(re) = this->GetVertexAttributes(re.first, &current_attribute, gibbs_args);
                        bottom_up(current_configuration);
                        energies.push_back(calculate_score(current_configuration));
                    }
                    std::vector<double> probs = neg_softmax(energies);
                    // if(this->GetVertexContent(re.first)->GetState().GetContent().compare(0, 3, "Fac") == 0 && a_pos == 3 && gr % 20 == 0)
                    // {   
                    //     std::cout << gr << std::endl;
                    //     std::cout << energies << std::endl;
                    //     std::cout << probs << std::endl;
                    // }
                    std::discrete_distribution<> dis(probs.begin(), probs.end());
                    unsigned index = dis(gen);
                    gibbs_args[1] = index;
                    current_configuration.at(re) = this->GetVertexAttributes(re.first, &current_attribute, gibbs_args);
                    bottom_up(current_configuration);
                }
            }
            configurations.push_back(current_configuration);
        }
        auto end_time = Clock::now();
        // if (gibbs_rounds != 0)
        //     std::cout << "Gibbs Sampling takes: " 
        //               << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() / gibbs_rounds
        //               << " milliseconds/round" << std::endl;
        return parse_tree;
    }

    template<class StateType, class AttributeType>    
    std::vector<Symbolic_State<StateType> > AOG<StateType, AttributeType>::GetTopLevelStates()
    {
        std::vector<Symbolic_State<StateType> > top_level_states;
        std::unordered_set<Symbolic_State<StateType> > sources;
        std::unordered_set<Symbolic_State<StateType> > results;
        for (Symbolic_Rule<StateType> rule : this->GetRules())
        {
            sources.insert(rule.GetSource());
            results.insert(rule.GetResults().begin(), rule.GetResults().end());
        }

        // check which source is not other sources' result
        for (Symbolic_State<StateType> source : sources)
        {
            if (results.find(source) == results.end())
                top_level_states.push_back(source);
        }
        return top_level_states;
    }

    template <class StateType, class AttributeType>
    double AOG<StateType, AttributeType>:: GrammarSize(double leaf_bias)
    {
        double grammar_size = 0;
        for(const auto& rule : this->GetRules())
        {
            ++grammar_size;
            for(const auto& state : rule.GetResults())
            {
                if(state.GetId() == -1)
                    grammar_size += leaf_bias;
                else    
                    ++grammar_size;
            }
        }
        return grammar_size;
    }

    template<class StateType, class AttributeType>
    void AOG<StateType, AttributeType>::TruncateGraph()
    {
        VertexId root_id = this->GetRoot();
        std::queue<VertexId> traverse_q;
        traverse_q.push(root_id);
        int count = 1;
        while (!traverse_q.empty())
        {
            VertexId cur_vertex_id = traverse_q.front();
            traverse_q.pop();
            std::vector<VertexId> children_id = this->ChildrenVertices(cur_vertex_id);
            for (VertexId child_id : children_id)
            {
                traverse_q.push(child_id);
            }

            if (children_id.size() == 1 && cur_vertex_id != this->GetRoot())
            {
                // find its parent, and connect child to parent
                std::vector<VertexId> parents_id = this->ParentsVertices(cur_vertex_id);
                // the reconnection is needed for all parents of cur_vertex_id
                for (VertexId p_id : parents_id)
                {
                    std::vector<Symbolic_State<StateType>> p_children_states;
                    // delete old edge parent to cur
                    std::vector<VertexId> p_children = this->ChildrenVertices(p_id);
                    int pos = -1;
                    std::vector<std::vector<std::shared_ptr<AOG_Edge> > > parent_to_p;
                    for (int i = 0; i < p_children.size(); i++)
                    {
                        if (p_children[i] == cur_vertex_id)
                            pos = i;
                        parent_to_p.push_back(this->GetEdgeContent(p_id, p_children[i]));
                        this->DeleteEdge(p_id, p_children[i]);
                        p_children_states.push_back(this->GetVertexContent(p_children[i])->GetState());
                    }

                    for (int i = 0; i < pos; i++)
                    {
                        std::shared_ptr<AOG_Edge> new_edge = std::make_shared<AOG_Edge>();                        
                        // Or node will not have multi-edge
                        if (!this->GetVertexContent(p_id)->IsNotOr())
                        {
                            if (parent_to_p[i].size() == 0)
                                throw std::exception();
                            if (parent_to_p[i].size() > 1)
                            {
                                std::cerr << "Multiedge\n";
                                throw std::exception();
                            }
                            new_edge->SetWeight(parent_to_p[i][0]->GetWeight());
                        }
                        // add new edge with weight if parent is Or
                        this->AddEdge(p_id, p_children[i], new_edge);
                    }

                    // repeat same operation for cur_vertex_id
                    // if parent is Or-node, copy its weight when creating the edge
                    // Or node will not have multi-edge
                    std::shared_ptr<AOG_Edge> new_edge = std::make_shared<AOG_Edge>();                        
                    if (!this->GetVertexContent(p_id)->IsNotOr())
                    {
                        if (parent_to_p[pos].size() == 0)
                            throw std::exception();
                        if (parent_to_p[pos].size() > 1)
                        {
                            std::cerr << "Multiedge\n";
                            throw std::exception();
                        }
                        new_edge->SetWeight(parent_to_p[pos][0]->GetWeight());
                    }
                    // add new edge with weight if parent is Or
                    this->AddEdge(p_id, children_id[0], new_edge);
                    for (int i = pos+1; i < p_children.size(); i++)
                    {
                        std::shared_ptr<AOG_Edge> new_edge = std::make_shared<AOG_Edge>();                                                
                        // Or node will not have multi-edge
                        if (!this->GetVertexContent(p_id)->IsNotOr())
                        {
                            if (parent_to_p[i].size() == 0)
                                throw std::exception();
                            if (parent_to_p[i].size() > 1)
                            {
                                std::cerr << "Multiedge\n";
                                throw std::exception();
                            }
                            new_edge->SetWeight(parent_to_p[i][0]->GetWeight());
                        }
                        // add new edge with weight if parent is Or
                        this->AddEdge(p_id, p_children[i], new_edge);
                    }
                    // Symbolic_Rule<StateType> to_delete(this->GetVertexContent(p_id)->GetState(), p_children_states);
                    // if(this->DeleteRule(to_delete))
                    // {
                    //     p_children_states[pos] = this->GetVertexContent(children_id[0])->GetState();
                    //     Symbolic_Rule<StateType> to_add(this->GetVertexContent(p_id)->GetState(), p_children_states);
                    //     this->AddRule(to_add);
                    // }
                }
                // delete old edge cur to child
                this->DeleteEdge(cur_vertex_id, children_id[0]);
                // if(this->GetVertexIdByState(this->GetVertexContent(cur_vertex_id)->GetState()) == cur_vertex_id)
                // {
                //     Symbolic_Rule<StateType> to_delete(this->GetVertexContent(cur_vertex_id)->GetState(),
                //                                        {this->GetVertexContent(children_id[0])->GetState()});
                //     this->DeleteRule(to_delete);
                // }
                // delete dummy node
                this->DeleteVertex(cur_vertex_id);
            }
        }
    }

    template<class StateType, class AttributeType>
    void AOG<StateType, AttributeType>::SimplifyGraph()
    {
        std::cout << "Simplifying..." << std::endl;

        //get max state id
        int stateId = 0;
        for(auto it : this->state_to_vertex_){
            if(it.first.GetId() > stateId){
                stateId = it.first.GetId();
            }
        }
        stateId++;
        
        auto delSameParAndChildState = [this](){
            std::queue<VertexId> q;
            q.push(this->GetRoot());

            bool OuterChanged = false;

            while(!q.empty()){
                VertexId curr = q.front();
                bool changed = false;
                unsigned childPos = 0;
                std::vector<VertexId> children = this->ChildrenVertices(curr);

                for(VertexId child : children){
                    
                    bool childIsNotOr = this->GetVertexContent(child)->IsNotOr();

                    if(!this->GetStateByVertexId(child).GetIsBasic() && 
                        childIsNotOr == this->GetVertexContent(curr)->IsNotOr() &&
                        this->ParentsVertices(child).size() == 1){
                        
                        OuterChanged = true;
                        changed = true;
                        this->DeleteEdge(curr, child);

                        if(childIsNotOr){
                            for(VertexId grandChild : this->ChildrenVertices(child)){
                                std::shared_ptr<AOG_Edge> tmpEdge = std::make_shared<AOG_Edge>();
                                this->AddEdge(curr, grandChild, tmpEdge, true, childPos);
                                childPos++;
                                if(!this->GetStateByVertexId(grandChild).GetIsBasic()){
                                    q.push(grandChild);
                                }
                            }
                            childPos--;
                        }
                        else {
                            for(VertexId grandChild : this->ChildrenVertices(child)){
                                auto it = std::find(children.begin(), children.end(), grandChild);
                                double weight = this->GetEdgeContent(child, grandChild)[0]->GetWeight();
                                if(it == children.end()){
                                    std::shared_ptr<AOG_Edge> tmpEdge = std::make_shared<AOG_Edge>(weight);
                                    this->AddEdge(curr, grandChild, tmpEdge, false); //or does not have multiedge
                                    if(!this->GetStateByVertexId(grandChild).GetIsBasic()){
                                        q.push(grandChild);
                                    }
                                }
                                else{ //rare case
                                    std::unordered_map<VertexId, double> m = this->GetOutEdgeWeights(curr, false);
                                    m[*it] += weight;
                                    this->SetOutEdgeWeights(curr, m);
                                }
                            }
                        }

                        this->DeleteVertex(child);
                    }

                    else if(!this->GetStateByVertexId(child).GetIsBasic()){
                        q.push(child);
                    }

                    childPos++;
                }
                
                if(!changed){
                    q.pop();
                }
            }

            return OuterChanged;
        };

        auto moveCommonChild = [this, &stateId]()
        {

            std::queue<VertexId> q;
            q.push(this->GetRoot());

            bool OuterChanged = false;

            while(!q.empty()){
                VertexId curr = q.front();
                q.pop();

                const std::vector<VertexId>& children = this->ChildrenVertices(curr);

                
                if(!this->GetVertexContent(curr)->IsNotOr() && children.size() > 1){
                    //find map from grandchildren in 1st position and last position to its parents (at child level wrt. curr)
                    std::unordered_map<VertexId, std::vector<VertexId>> fstGrandChildren, sndGrandChildren;
                    for (VertexId child : children){
                        const std::vector<VertexId>& grandChildren = this->ChildrenVertices(child);
                        if(grandChildren.empty()){
                            continue;
                        }
                        VertexId fstGrandChild = grandChildren.front();
                        VertexId sndGrandChild = grandChildren.back();
                        auto fst_it = fstGrandChildren.find(fstGrandChild);
                        if(fst_it == fstGrandChildren.end()){
                            fstGrandChildren[fstGrandChild] = {child};
                        }
                        else{
                            fst_it->second.push_back(child);
                        }
                        auto snd_it = sndGrandChildren.find(sndGrandChild);
                        if(snd_it == sndGrandChildren.end()){
                            sndGrandChildren[sndGrandChild] = {child};
                        }
                        else{
                            snd_it->second.push_back(child);
                        }
                    }
                    
                    //find map from child level to their common grandchildren
                    std::map<std::vector<VertexId>, std::pair<VertexId, VertexId>> Children2Grand;
                    for(const auto& x : fstGrandChildren){
                        if(x.second.size() >= 2)
                            Children2Grand[x.second] = {x.first, 0};
                    }
                    for(const auto& x : sndGrandChildren){
                        if(x.second.size() >= 2){
                            auto it = Children2Grand.find(x.second);
                            if(it == Children2Grand.end()){
                                Children2Grand[x.second] = {0, x.first};
                            }
                            else if(it->second.first != x.first){
                                it->second.second = x.first;
                            }
                        }
                    }
                    
                    if(Children2Grand.empty()){
                        for (VertexId child : children){
                            if(!this->GetStateByVertexId(child).GetIsBasic()){
                                q.push(child);
                            }
                        }
                        continue;
                    }

                    //if conflict in child level occurs, keep the child level with larger size
                    std::vector<std::pair<std::vector<VertexId>, std::pair<VertexId, VertexId>>> ch2g; 
                    std::copy(Children2Grand.begin(), Children2Grand.end(), back_inserter(ch2g));
                    std::sort(ch2g.begin(), ch2g.end(), 
                            [](const std::pair<std::vector<VertexId>, std::pair<VertexId, VertexId>> &a, const std::pair<std::vector<VertexId>, std::pair<VertexId, VertexId>> & b)
                            {return a.first.size() < b.first.size();});
                    
                    std::unordered_set<VertexId> existVertex;
                    for(auto it = ch2g.begin(); it != ch2g.end();){
                        bool found = false;
                        for(VertexId x : it->first){
                            if(existVertex.find(x) != existVertex.end()){
                                it = ch2g.erase(it);
                                found = true;
                                break;
                            }
                        }
                        if(found){
                            continue;
                        }
                        for(VertexId x : it->first){
                            existVertex.insert(x);
                        }
                        ++it;
                    }

                    //push no common grandchildren child to the queue
                    for(VertexId child : children){
                        if(existVertex.find(child) == existVertex.end() && !this->GetStateByVertexId(child).GetIsBasic()){
                            q.push(child);
                        }
                    }
                    
                    for(const auto& x : ch2g){
                        //test for children, should have no parents other than curr
                        const std::vector<VertexId>& andVertices = x.first;
                        VertexId fstGrandChild = x.second.first;
                        VertexId sndGrandChild = x.second.second;

                        bool noOtherParent = true;
                        for(VertexId child : andVertices){
                            std::vector<VertexId> otherParents = this->ParentsVertices(child);
                            if(otherParents.size() != 1){
                                noOtherParent = false;
                                break;
                            }
                        }

                        if(!noOtherParent){
                            for (VertexId child : andVertices){
                                if(!this->GetStateByVertexId(child).GetIsBasic()){
                                    q.push(child);
                                }
                            }
                            continue;
                        }

                        // std::cout << "[simplify] pass all tests" << std::endl;

                        //separate cases for:
                        //1. all and nodes share common children, or
                        //2. some of the and nodes do
                        //the separation of cases is for the purpose of keep state id unchanged after truncate() for case 1
                        if(andVertices.size() == children.size()){ //case 1
                            std::shared_ptr<Symbolic_State<StateType> > dummy_state = std::make_shared<Symbolic_State<StateType>>(stateId++);
                            std::shared_ptr<AOG_Vertex<StateType, AttributeType> > dummy = std::make_shared<AOG_Vertex<StateType, AttributeType>>(*dummy_state, true, false);
                            unsigned dummy_id = this->AddVertex(dummy);

                            //connect dummy to parents
                            for(VertexId parent : this->ParentsVertices(curr)){
                                if(!this->GetVertexContent(parent)->IsNotOr()){ //need to consider weight, no multiedge case, and no position requirement
                                    std::shared_ptr<AOG_Edge> tmpEdge = std::make_shared<AOG_Edge>(this->GetEdgeContent(parent, curr)[0]->GetWeight());
                                    this->AddEdge(parent, dummy_id, tmpEdge, false);
                                }
                                else{
                                    std::vector<VertexId> childrenOfPar = this->ChildrenVertices(parent);
                                    std::vector<unsigned> pos;
                                    auto iter = childrenOfPar.begin();
                                    auto iter_begin = childrenOfPar.begin();
                                    while ((iter = std::find(iter, childrenOfPar.end(), curr)) != childrenOfPar.end()){
                                        pos.push_back(iter - iter_begin);
                                        iter++;
                                    }
                                    this->DeleteEdge(parent, curr);
                                    for(unsigned innerPos : pos){
                                        std::shared_ptr<AOG_Edge> tmpEdge = std::make_shared<AOG_Edge>();
                                        this->AddEdge(parent, dummy_id, tmpEdge, true, innerPos);
                                    }
                                }
                            }

                            if(curr == this->GetRoot()){
                                this->root_id_ = dummy_id;
                                this->GetVertexContent(curr)->SetIsRoot(false);
                                this->GetVertexContent(dummy_id)->SetIsRoot(true);
                            }

                            //connect dummy to fstGrandChild (if shared), curr, and lstGrandChild (if shared)
                            if(fstGrandChild != 0){
                                for(VertexId child : children){
                                    this->DeleteEdge(child, fstGrandChild);
                                }
                                std::shared_ptr<AOG_Edge> tmpEdge = std::make_shared<AOG_Edge>();
                                this->AddEdge(dummy_id, fstGrandChild, tmpEdge, true);
                                if(!this->GetStateByVertexId(fstGrandChild).GetIsBasic()){
                                    q.push(fstGrandChild);
                                }
                            }

                            std::shared_ptr<AOG_Edge> tmpEdge = std::make_shared<AOG_Edge>();
                            this->AddEdge(dummy_id, curr, tmpEdge, true);
                            q.push(curr);

                            if(sndGrandChild != 0){
                                for(VertexId child : children){
                                    this->DeleteEdge(child, sndGrandChild);
                                }
                                std::shared_ptr<AOG_Edge> tmpEdge = std::make_shared<AOG_Edge>();
                                this->AddEdge(dummy_id, sndGrandChild, tmpEdge, true);
                                if(!this->GetStateByVertexId(sndGrandChild).GetIsBasic()){
                                    q.push(sndGrandChild);
                                }
                            }

                            OuterChanged = true;
                        }
                        
                        else{ //case 2
                            std::shared_ptr<Symbolic_State<StateType> > dummy_and_state = std::make_shared<Symbolic_State<StateType>>(stateId++);
                            std::shared_ptr<AOG_Vertex<StateType, AttributeType> > dummy_and = std::make_shared<AOG_Vertex<StateType, AttributeType>>(*dummy_and_state, true, false);
                            unsigned dummy_and_id = this->AddVertex(dummy_and);

                            std::shared_ptr<Symbolic_State<StateType> > dummy_or_state = std::make_shared<Symbolic_State<StateType>>(stateId++);
                            std::shared_ptr<AOG_Vertex<StateType, AttributeType> > dummy_or = std::make_shared<AOG_Vertex<StateType, AttributeType>>(*dummy_or_state, false, false);
                            unsigned dummy_or_id = this->AddVertex(dummy_or);

                            //connect dummy_and to curr, disconnect selected children to curr, connect children to dummy_or
                            double dummy_and_weight = 0;
                            for(VertexId andVtx : andVertices){
                                double weight = this->GetEdgeContent(curr, andVtx)[0]->GetWeight();
                                this->DeleteEdge(curr, andVtx);
                                std::shared_ptr<AOG_Edge> dummy_or_edge = std::make_shared<AOG_Edge>(weight);
                                this->AddEdge(dummy_or_id, andVtx, dummy_or_edge, false);
                                dummy_and_weight += weight;
                            }
                            std::shared_ptr<AOG_Edge> dummy_and_edge = std::make_shared<AOG_Edge>(dummy_and_weight);
                            this->AddEdge(curr, dummy_and_id, dummy_and_edge, false);

                            //connect dummy to fstGrandChild (if shared), curr, and lstGrandChild (if shared)
                            if(fstGrandChild != 0){
                                for(VertexId andVtx : andVertices){
                                    this->DeleteEdge(andVtx, fstGrandChild);
                                }
                                std::shared_ptr<AOG_Edge> tmpEdge = std::make_shared<AOG_Edge>();
                                this->AddEdge(dummy_and_id, fstGrandChild, tmpEdge, true);
                                if(!this->GetStateByVertexId(fstGrandChild).GetIsBasic()){
                                    q.push(fstGrandChild);
                                }
                            }

                            std::shared_ptr<AOG_Edge> dummy_and_to_or_edge = std::make_shared<AOG_Edge>();
                            this->AddEdge(dummy_and_id, dummy_or_id, dummy_and_to_or_edge, true);
                            q.push(dummy_or_id);

                            if(sndGrandChild != 0){
                                for(VertexId andVtx : andVertices){
                                    this->DeleteEdge(andVtx, sndGrandChild);
                                }
                                std::shared_ptr<AOG_Edge> tmpEdge = std::make_shared<AOG_Edge>();
                                this->AddEdge(dummy_and_id, sndGrandChild, tmpEdge, true);
                                if(!this->GetStateByVertexId(sndGrandChild).GetIsBasic()){
                                    q.push(sndGrandChild);
                                }
                            }

                            OuterChanged = true;
                        }
                    }
                }
                else{
                    for (VertexId child : children){
                        if(!this->GetStateByVertexId(child).GetIsBasic()){
                            q.push(child);
                        }
                    }
                }
            }

            return OuterChanged;
        };

        while(true){
            bool changed = delSameParAndChildState();
            bool changed2 = moveCommonChild();
            changed = changed || changed2;
            if(!changed){
                break;
            }
            this->TruncateGraph();
        }
    }

    template<class StateType, class AttributeType>
    bool AOG<StateType, AttributeType>::Merge(const AOG<StateType, AttributeType> &targetAOG) 
    {
        auto newRules = targetAOG.GetRules();

        for(auto rule : newRules){
            this->AddRule(rule);
        }

        return true;
    }

    template<class StateType, class AttributeType>
    void AOG<StateType, AttributeType>::DeleteNoParentRules(const Symbolic_State<StateType>& src_state)
    {
        VertexId vid = this -> GetVertexIdByState(src_state);
        
        //if the passed in vertex has parent, is root, or the state is a leaf, directly return
        if(this->ParentsVertices(vid).size() || src_state.GetIsBasic() || vid == this->root_id_)
            return;
        std::vector<VertexId> children_vtx = this->ChildrenVertices(vid);
        // if it is an and-node
        if(this->GetVertexContent(vid)->IsNotOr())
        {
            //delete curent rule
            std::vector<Symbolic_State<StateType> > children_states;
            for(auto vertex : children_vtx )
                children_states.push_back(this->GetStateByVertexId(vertex));
            
            if(!this->DeleteRule(Symbolic_Rule<StateType>(src_state ,children_states)))
                std::cerr<<"fail to delete rules in and-node:\n";
            
            //recursively delete child states if they have no parents
            for(auto child_state : children_states)
                this->DeleteNoParentRules(child_state);
        }
        //if it is an or-node
        else
        {
            std::vector< std::vector< Symbolic_State<StateType> > > grand_children_states;
            //first create all rules
            for(auto child : children_vtx)            
            {
                std::vector<VertexId> grand_children_vtx = this->ChildrenVertices(child);
                std::vector<Symbolic_State<StateType> > seq;
                for(auto grand_child : grand_children_vtx )
                    seq.push_back(this->GetStateByVertexId(grand_child));
                grand_children_states.push_back(seq);
            }
            //delete all rules
            for(int i = 0; i < grand_children_states.size(); ++i)
            {
                if(!this->DeleteRule(Symbolic_Rule<StateType>(src_state ,grand_children_states[i])))
                    std::cerr<<"fail to delete rules in or-node\n";
                for(auto state : grand_children_states[i])
                    this->DeleteNoParentRules(state);
            }   
        }
    }

    template<class StateType, class AttributeType>
    void AOG<StateType, AttributeType>::Visualize(std::string dir, std::string filename, bool truncate)
    {
        std::cout << "PRINTED GRAPH\n";
        std::shared_ptr<AOG<StateType, AttributeType> > printed_graph =
            std::make_shared<AOG<StateType, AttributeType> >(*this);

        if (truncate)
            printed_graph->TruncateGraph();

        std::string dir_copy = dir;
        std::string filename_copy = filename;
        std::ofstream file;
        if (dir_copy.empty())
            dir_copy = filename_copy;
        else
            dir_copy = dir_copy.append("/").append(filename_copy);

        file.open(dir_copy, std::ofstream::out|std::ofstream::trunc);
        if (file.is_open()) {

            VertexId root_id = printed_graph->GetRoot();
            std::queue<VertexId> q;
            q.push(root_id);
            std::queue<unsigned> level_q;
            level_q.push(1);

            std::unordered_set<VertexId> plot_visited;

            while (!q.empty())
            {
                VertexId cur_vertex = q.front();
                unsigned level = level_q.front();
                q.pop();
                level_q.pop();

                std::vector<VertexId> children_vertices_id = printed_graph->ChildrenVertices(cur_vertex);

                if(!printed_graph->GetVertexContent(cur_vertex))
                {
                    int a;
                }
                bool IsNotOr = printed_graph->GetVertexContent(cur_vertex)->IsNotOr();
                std::unordered_map<VertexId, double> weights;
                Symbolic_State<StateType> cur_vertex_state = printed_graph->GetStateByVertexId(cur_vertex);
                
                if (!IsNotOr)
                    weights = printed_graph->GetOutEdgeWeights(cur_vertex, true);
                    
                for (int i = 0; i < children_vertices_id.size(); i++)
                {
                    Symbolic_State<StateType> child_state = printed_graph->GetStateByVertexId(children_vertices_id[i]);
                    if (IsNotOr)
                    {    
                        file << cur_vertex << "," << children_vertices_id[i]
                             << "," << cur_vertex_state.GetContent() << "_"
                             << cur_vertex_state.GetId() << "," << child_state.GetContent()
                             << "_" << child_state.GetId() << ",," << i+1 << "," << level << "\n";
                    }
                    else
                    {
                        file << cur_vertex << "," << children_vertices_id[i]
                             << "," << cur_vertex_state.GetContent() << "_"
                             << cur_vertex_state.GetId() << "," << child_state.GetContent()
                             << "_" << child_state.GetId() << "," << weights[children_vertices_id[i]]
                             << "," << i+1 <<"," << level << "\n";
                    }
                    q.push(children_vertices_id[i]);
                    level_q.push(level+1);
                }
            }
            file.close();
        }
        else
            std::cout << "Unable to open " << dir_copy << std::endl;
    }

    template<class StateType, class AttributeType>
    void AOG<StateType, AttributeType>::SaveGraph(std::string path, std::string name)
    {
        std::string file_name;

        file_name = path + '/' + name;
        
        std::ofstream file;
        file.open(file_name, std::ofstream::out|std::ofstream::trunc);

        if (file.is_open())
        {
            std::cerr << "File opened. \n";
            for (const Symbolic_Rule<StateType>& rule : this->all_rules_)
            {
                // size of result, souce id, source content, [child id, child content]*
                file << rule.GetResults().size();
                std::vector<Symbolic_State<StateType> > states = rule.GetResults();
                VertexId source_vtx_id = this->GetVertexIdByState(rule.GetSource());
                bool IsNotOr = this->GetVertexContent(source_vtx_id)->IsNotOr();
                if(!IsNotOr){
                    std::vector<VertexId> children_vector;
                    for(const Symbolic_State<StateType>& state : states)
                        children_vector.push_back(this->GetVertexIdByState(state));
                    VertexId tmpDummy;
                    bool found = false;
                    for (VertexId dummy : this ->ChildrenVertices(source_vtx_id))
                    {
                        if (this->ChildrenVertices(dummy) == children_vector)
                        {
                            found = true;
                            tmpDummy = dummy;
                            break;
                        }
                    }
                    if (!found)
                    {
                        std::cerr << "rule not found" << std::endl;
                        throw std::exception();
                    }
                    double weight = this->GetOutEdgeWeights(source_vtx_id, false)[tmpDummy];
                    if (weight > 0)
                        file << "," << weight;
                    else{
                        std::cerr << "0 or negative weight is found !!" << std::endl;
                        throw std::exception();
                    }
                }
                else{
                    file << "," << 0;
                }
                file << "," << rule.GetSource().GetId() << "," << rule.GetSource().GetContent();
                for (int i = 0; i < states.size(); i++)
                {
                    file << "," << states[i].GetId() << "," << states[i].GetContent();
                }
                file << "\n";
            }

            file.close();
        }
        else
            std::cout << "Unable to open " << "learned_tree.txt" << std::endl;
    }
}

#endif//AOG_LIB_AOG_H
