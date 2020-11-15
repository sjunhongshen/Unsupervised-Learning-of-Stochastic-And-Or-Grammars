//
// Created by Luyao Yuan on 18/2/22.
//

#ifndef AOG_LIB_CONTEXT_MATRIX_H
#define AOG_LIB_CONTEXT_MATRIX_H

#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/multi_index/member.hpp>

#include "AOG.h"

template <class T>
using SequenceType = std::vector<AOG_LIB::Symbolic_State<T> >;

template <class T>
using ContextType = std::vector<SequenceType<T> >;

template<class T>
bool operator==(const SequenceType<T>& rhs, const SequenceType<T>& lhs)
{
    if(rhs.size() != lhs.size()) return false;
    unsigned i = 0;
    while(i < rhs.size())
        if(!(rhs[i] == lhs[i++]))
            return false;
    return true;
}

template <class T>
bool operator==(const ContextType<T>& rhs, const ContextType<T>& lhs)
{
    if(rhs.size() != lhs.size())
        return false;
    
    if(!rhs.size() || !lhs.size() )
    {
        std::cerr<<"one of the contexts is empty!\n";
        throw std::exception();
    }

    unsigned i = 0;
    while(i < rhs.size())
        if(!(rhs[i] == lhs[i++]))
            return false;
    return true;

}

template <class StateType>
size_t hash_value(const SequenceType<StateType> &seq)
{
    // using std::string;
    size_t vec_seed = boost::hash_range(seq.begin(),seq.end());
    return vec_seed;
}

namespace std
{
    template<class StateType>
    struct hash<ContextType<StateType> >
    {
        size_t operator()(const ContextType<StateType> &context) const noexcept
        {

            using std::string;
            size_t vec_seed = boost::hash_range(context.begin(), context.end());
            // size_t vec_seed_2 = boost::hash_range(context.second.begin(), context.second.end());

            // Compute individual hash values for first,
            // second and third and combine them using XOR
            // and bit shifting:

            return vec_seed;
        }
    };
}

namespace AOG_LIB
{
    //An entry in the context matrix, the data can be recovered by
    //[context_[0], configuration[0], context_[1], ...]
    template <class StateType>
    struct CM_Entry
    {
        ContextType<StateType> context_;
        std::vector<SequenceType<StateType> >  configuration_;
        unsigned count_;
    };
}

template <class T>
using Context_Matrix = boost::multi_index::multi_index_container
<
    AOG_LIB::CM_Entry<T>,
    boost::multi_index::indexed_by
    <
        boost::multi_index::hashed_non_unique //query by context
        <
            boost::multi_index::member<AOG_LIB::CM_Entry<T>, ContextType<T>, &AOG_LIB::CM_Entry<T>::context_>
        >,
        boost::multi_index::hashed_non_unique //query by configuration
        <
            boost::multi_index::member<AOG_LIB::CM_Entry<T>, std::vector<SequenceType <T> >, &AOG_LIB::CM_Entry<T>::configuration_>
        >,
	    boost::multi_index::hashed_non_unique //query by count, may not be necessary, remove if not
	    <
		    boost::multi_index::member<AOG_LIB::CM_Entry<T>, unsigned , &AOG_LIB::CM_Entry<T>::count_>
        >,
        boost::multi_index::random_access<> //enable vector like access, .size() can get memory usage
    >
>;

#endif //AOG_LIB_CONTEXT_MATRIX_H