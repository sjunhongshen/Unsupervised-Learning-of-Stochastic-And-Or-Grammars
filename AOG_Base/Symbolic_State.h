//
// Created by Luyao Yuan on 17/10/23.
//

#ifndef AOG_LIB_SYMBOLIC_STATE_H
#define AOG_LIB_SYMBOLIC_STATE_H

#include <functional>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/functional/hash.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

namespace AOG_LIB
{
    /*
        * This class defines symbolic states, which are the semantic meanings of each AOG vertex in AOG graph.
        */
    template <class StateType>
    class Symbolic_State
    {
        friend class boost::serialization::access;
        // When the class Archive corresponds to an output archive, the
        // & operator is defined similar to <<.  Likewise, when the class Archive
        // is a type of input archive the & operator is defined similar to >>.
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & state_content_;
            ar & is_basic_;
            ar & id_;
        }
        StateType state_content_; // the content of the symbolic state
        bool is_basic_;           // true if the symbolic state is used in a leaf node
        int id_;                  //if the state has semantic meaning, then -1, else other integers used to distinguish states with no semantic meaning

    public:
        Symbolic_State();
        /**
         * constructor for state with no actual semantic meaning
         */
        Symbolic_State(int);
       
        /**
         * constructs a symbolic state with a given state and nodes' properties
         * id_ is -1 in this case
         * 
         * @param content: the semantic meaning of the state
         * @param is_leaf: true if the symbolic state is used in a leaf node
         */
        Symbolic_State(const StateType&, bool);

        /**
         * Check is the state is a basic state 
         * 
         * @returns true if the symbolic state is used in a leaf node.
         */
        bool GetIsBasic() const;

        /**
         * Get the content of the state
         * 
         * @returns the content(the semantic meaning) of the symbolic state.
         */
        const StateType &GetContent() const;

        /**
         * Get the id of the state
         * 
         * @returns the id of the symbolic state.
         */
        int GetId() const { return this->id_; };

        /**
         * This function compares two symbolic states
         * @param rhs: the other symbolic state to be compared with
         * @returns true if the content, is_basic of the two symbolic states are all the same
         */
        bool operator==(const Symbolic_State<StateType> &) const;
        bool operator!=(const Symbolic_State<StateType> &) const;
    };
    template <class StateType>
    size_t hash_value(const AOG_LIB::Symbolic_State<StateType> &state)
    {

        using std::string;

        // Compute individual hash values for first,
        // second  and combine them using XOR
        // and bit shifting:
        size_t seed = 0;
        boost::hash_combine(seed, state.GetContent());
        boost::hash_combine(seed, state.GetIsBasic());
        boost::hash_combine(seed, state.GetId());
        return seed;
    }

    template <class StateType>
    size_t hash_value(const std::vector<AOG_LIB::Symbolic_State<StateType> > &seq)
    {
        // using std::string;
        size_t vec_seed = boost::hash_range(seq.begin(),seq.end());
        return vec_seed;
    }
}

namespace std
{
    template <class StateType>
    struct hash<AOG_LIB::Symbolic_State<StateType> >
    {
        size_t operator()(const AOG_LIB::Symbolic_State<StateType> &state) const noexcept
        {
            boost::hash<AOG_LIB::Symbolic_State<StateType> > hasher;
            return hasher(state);
        }
    };

    template<class StateType>
    struct hash<vector<AOG_LIB::Symbolic_State<StateType> > >
    {
        size_t operator()(const vector<AOG_LIB::Symbolic_State<StateType> > &seq) const noexcept
        {

            // using std::string;
            boost::hash<vector<AOG_LIB::Symbolic_State<StateType> > > seq_hasher;
            return seq_hasher(seq);
        }
    };
}

template <class StateType>
AOG_LIB::Symbolic_State<StateType>::Symbolic_State()
    : id_(-1), is_basic_(true)
{}

//StateType has to support default constructor to use this function
//create a state with no semantic meaning
template <class StateType>
AOG_LIB::Symbolic_State<StateType>::Symbolic_State(int id)
    : id_(id), is_basic_(false)
{}

template <class StateType>
AOG_LIB::Symbolic_State<StateType>::Symbolic_State(const StateType &content, bool is_leaf)
    : state_content_(content), is_basic_(is_leaf), id_(-1)
{}

template <class StateType>
bool AOG_LIB::Symbolic_State<StateType>::GetIsBasic() const
{return this->is_basic_;}

template <class StateType>
const StateType &AOG_LIB::Symbolic_State<StateType>::GetContent() const
{return this->state_content_;}

template <class StateType>
bool AOG_LIB::Symbolic_State<StateType>::operator==(const Symbolic_State<StateType> &rhs) const
{
    return (this->is_basic_ == rhs.GetIsBasic() &&
            this->state_content_ == rhs.GetContent() &&
            this->id_ == rhs.GetId());
}

template <class StateType>
bool AOG_LIB::Symbolic_State<StateType>::operator!=(const Symbolic_State<StateType> &rhs) const
{
    return !(*this == rhs);
}

template <class StateType>
std::vector<std::vector<AOG_LIB::Symbolic_State<StateType> > > FileParser(std::string filename)
{
    std::ifstream file(filename);
	if(!file.good())
	{
		std::cout << "Cannot open file " << filename << std::endl;
		throw std::exception();
	}
    std::string line;
    std::vector<std::vector<AOG_LIB::Symbolic_State<StateType> > > ParsedVector;

    while (std::getline(file, line))
    {
        std::istringstream iss(line);

        std::vector<AOG_LIB::Symbolic_State<StateType> > subVector;
        do
        {
            std::string subs;
            iss >> subs;        
            AOG_LIB::Symbolic_State<StateType> subState(subs, true);
            subVector.push_back(subState);
        } while (iss);
        subVector.pop_back();
        ParsedVector.push_back(subVector);
    }
    return ParsedVector;
}

#endif //AOG_LIB_SYMBOLIC_STATE_H
