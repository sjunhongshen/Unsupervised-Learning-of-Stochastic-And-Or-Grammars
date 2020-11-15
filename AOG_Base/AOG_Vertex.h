//
// Created by yuanluyao on 10/23/17.
//

#ifndef AOG_LIB_AOG_VERTEX_H
#define AOG_LIB_AOG_VERTEX_H

#include <string>
#include <functional>
#include "Symbolic_State.h"
#include "../Core/Graph.hpp"
namespace AOG_LIB
{
    //forward declaration
    template <class StateType, class AttributeType> class AOG;
    class AOG_Edge;

    /* This class defines vertices in an AOG graph.
     * Mutators in this class are set as private to prevent unwanted changes of vertices' properties from users
     * This class declares T_AOG as friend class to enable T_AOG to change its vertices' properties.
     */
    template<class StateType, class AttributeType>
    class AOG_Vertex
    {
        bool is_not_or_; // true if the AOG Vertex is an And-node or containing a basic state
        bool is_root_;// true if the AOG Vertex is a root node
        Symbolic_State<StateType> state_; // The Symbolic state contained in the AOG Vertex
        //attribute update function
        std::function<std::vector<AttributeType>(const AOG_Vertex&,
            const std::vector<std::vector<AttributeType> >*, void*)> attribute_func_ = 0;
        //attributes of this vertex
        std::vector<std::vector<AttributeType> > attributes_ranges_;
        std::function<std::vector<std::vector<AttributeType> >(const AOG_Vertex&, void*)> attributes_ranges_func_ = 0;
        void* attributes_ranges_param_;


        //unary potential function
        std::function<double(const std::vector<AttributeType>&, void*)> potential_func_ = 0;
        void* potential_params_ = 0;

        //gradient function for unary potential function
        std::function<void(const std::vector<AttributeType>&, void*, std::vector<double>&)> potential_grad_func_ = 0;

        friend class AOG<StateType, AttributeType>;

        /* This function sets the symbolic state of a AOG Vertex.
         * @param:
         *      state: The symbolic state to be set
         * @return: void
         */
        void SetState(const Symbolic_State<StateType> &state){this->state_ = state;};
        
        /* This function sets a vertex to be an and-node or an or-node in AOG graph
         * @param:
         *      is_and: a boolean value that specifies whether the vertex is an and-node
         * @return: void
         */
        void SetIsNotOr(const bool is_not_or){this->is_not_or_ = is_not_or;};

        /**
         * This function sets a vertex's attribute calculate function
         * @param:
         *      func: the function to be used as the attribute function
         * @return:
         *      void
        */
        void SetAttributeFunc(const std::function<std::vector<AttributeType>(const AOG_Vertex&,
                                                         const std::vector<std::vector<AttributeType> >*, void*)> func)
        {this->attribute_func_ = func;};

        /**
         * This function sets a vertex's attribute ranges calculate function
         * @param:
         *      func: the function to be used as the attribute ranges function
         *      param: the param used for this function
         * @return:
         *      void
        */
        void SetAttributesRangesFunc
        (std::function<std::vector<std::vector<AttributeType> >(const AOG_Vertex&, void*)> func, void* param)
        {
            this->attributes_ranges_func_ = func;
            this->attributes_ranges_param_ = param;
        };

        /**
         * This function sets a vertex's potential function
         * @param:
         *      func: the function to be used as the potential function
         * @return:
         *      void
        */
        void SetPotentialFunc(const std::function<double(const std::vector<AttributeType>&, void*)> func)
        {this->potential_func_ = func;};

        /**
         * This function sets a vertex's potential gradient calculate function
         * @param:
         *      func: the function to be used as the gradient function of the potential function
         * @return:
         *      void
        */
        void SetPotentialGradFunc(const std::function<void(const std::vector<AttributeType>&, void*, std::vector<double>&)> func)
        {this->potential_grad_func_ = func;};

        /**
         * This function sets a vertex to be a root or not root
         * @param:
         *      is_root: a boolean value that specifies whether the vertex is the root
         * @return:
         *      void
        */
        void SetIsRoot(bool is_root){this->is_root_ = is_root;};

        /**
         * Get current attribute of this vertex
         * 
         * @param neighbors_attrs: neighbor attributes provided by the AOG
         * @param args: any relevant variables
         * @returns current attributes of this vertex
         */
        std::vector<AttributeType> GetAttributes(const std::vector<std::vector<AttributeType> >* neighbors_attrs, void *args) const;

        /**
         * Set attributes ranges for this vertex, every vector is a range for an attribute
         * 
         * @param: the range matrix
         */
        void SetAttributesRanges(const std::vector<std::vector<AttributeType>> &);

    public:
        /**
         * Construct an AOG Vertex from a symbolic state along with the node's properties.
         * @param state: the symbolic state used to construct the vertex
         * @param bool is_and: specify whether the vertex is an and-node or an or-node
         * @param bool is_root: specify whether the vertex is a root node
         */
        AOG_Vertex(const Symbolic_State<StateType> &, bool, bool,
                   const std::function<std::vector<AttributeType>(const AOG_Vertex&,
                                              const std::vector<std::vector<AttributeType> >*, void*)> a_func = 0);

        /**
         * Get the current symbolic state of the AOG vertex
         * 
         * @returns the current symbolic state of the AOG vertex
         */
        const Symbolic_State<StateType>& GetState() const;

        /**
         * Get the attribute's range of the AOG vertex to sample from
         * 
         * @returns the attribute's range of the AOG vertex
         */
        std::vector<std::vector<AttributeType>> GetAttributesRange() const;

        /**
         * Check if current node is the AND node
         * 
         * @returns a boolean value indicating whether the vertex is an and-node.
         */
            bool IsNotOr() const;

        /**
         * Check if current node is the root node
         * 
         * @returns a boolean value indicating whether the vertex is a root node.
         */
        bool IsRoot() const;

        /**
         * Get current attribute function of this vertex
         * 
         * @returns current attribute function of this vertex
         */
        std::function<std::vector<AttributeType>(const AOG_Vertex &, const std::vector<std::vector<AttributeType> >*, void*)>
        GetAttributeFunc() const;

        /**
         * Get potential function of this vertex
         * 
         * @returns current potential function of this vertex
         */
        std::function<double(const std::vector<AttributeType> &, void*)>
        GetPotentialFunc() const;

        /**
         * Get gradient function of the potential function of this vertex
         * 
         * @returns current gradient function of this vertex
         */
        std::function<void(const std::vector<AttributeType> &, void*, std::vector<double>&)>
        GetPotentialGradFunc() const;

        /**
         * This function sets a vertex's attribute calculate function
         * @param:
         *      param: the parameter of the potential function to be set
         * @return:
         *      void
        */
        void SetPotentialParam(void* param)
        {this->potential_params_ = param;};

        /**
         * This function sets a vertex's attribute calculate function
         * @return:
         *      the parameter of the potential function
        */
        void* GetPotentialParam() const
        {return this->potential_params_;};

        /**
         * This function updates this vertex's attribute ranges
         * @return:
         *      the parameter of the potential function
        */
        void UpdateAttributesRanges()
        {
            if(this->attributes_ranges_func_)
                this->attributes_ranges_ = this->attributes_ranges_func_(*this, this->attributes_ranges_param_);
        };
    };


    template<class StateType, class AttributeType>
    AOG_Vertex<StateType, AttributeType>::AOG_Vertex(const Symbolic_State<StateType> &state, bool is_not_or, bool is_root,
                                                     const std::function<std::vector<AttributeType>(const AOG_Vertex&,
                                                        const std::vector<std::vector<AttributeType> >*, void*)> a_func)
    :state_(state), is_not_or_(is_not_or), is_root_(is_root), attribute_func_(a_func)
    {}

    template<class StateType, class AttributeType>
    std::vector<AttributeType> AOG_Vertex<StateType, AttributeType>::
    GetAttributes(const std::vector<std::vector<AttributeType> >* neighbors_attrs, void* args) const
    {
        if(!this->attribute_func_)
        {
            std::cerr << "Vertex doesn't have attribute function\n";
            exit(1);
        }
        return this->attribute_func_(*this, neighbors_attrs, args);
    }

    template <class StateType, class AttributeType>
    void AOG_Vertex<StateType, AttributeType>::SetAttributesRanges(const std::vector<std::vector<AttributeType>> & attrs_ranges)
    {
        this->attributes_ranges_ = attrs_ranges;
    }

    template<class StateType, class AttributeType>
    const Symbolic_State<StateType>& AOG_Vertex<StateType, AttributeType>::GetState() const
    { return this->state_; }

    template<class StateType, class AttributeType>
    bool AOG_Vertex<StateType, AttributeType>::IsNotOr() const
    { return this->is_not_or_; }

    template<class StateType, class AttributeType>
    bool AOG_Vertex<StateType, AttributeType>::IsRoot() const
    { return this->is_root_; }

    template<class StateType, class AttributeType>
    std::vector<std::vector<AttributeType> > AOG_Vertex<StateType, AttributeType>::GetAttributesRange() const
    { return this->attributes_ranges_;}

    template <class StateType, class AttributeType>
    std::function<std::vector<AttributeType>(const AOG_Vertex<StateType, AttributeType> &,
                  const std::vector<std::vector<AttributeType> >*, void*)>
    AOG_Vertex<StateType, AttributeType>::GetAttributeFunc() const
    {return this->attribute_func_;}

    template <class StateType, class AttributeType>
    std::function<double(const std::vector<AttributeType>&, void*)>
    AOG_Vertex<StateType, AttributeType>::GetPotentialFunc() const
    {return this->potential_func_;}

    template <class StateType, class AttributeType>
    std::function<void(const std::vector<AttributeType>&, void*, std::vector<double>&)>
    AOG_Vertex<StateType, AttributeType>::GetPotentialGradFunc() const
    {return this->potential_grad_func_;}
}
#endif //AOG_LIB_AOG_VERTEX_H
