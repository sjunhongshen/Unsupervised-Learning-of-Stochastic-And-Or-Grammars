#ifndef AOG_LIB_EARLEY_EVALUATION_H
#define AOG_LIB_EARLEY_EVALUATION_H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <iterator>
#include <list>
#include <sstream>
#include <set>
#include <algorithm>
#include <utility>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <functional>
#include <boost/functional/hash.hpp>

#include "Symbolic_Rule.h"
#include "Symbolic_State.h"
#include "AOG.h"


#ifndef CONTAINER_HASH
#define CONTAINER_HASH
template <typename Container> // (hasher for containers) we make this generic for any container
struct container_hash {
    std::size_t operator()(Container const& c) const {
        return boost::hash_range(c.begin(), c.end());
    }
};
#endif

namespace AOG_LIB_UTIL
{

    using AOG_LIB::Symbolic_State;
    using AOG_LIB::Symbolic_Rule;
    using AOG_LIB::AOG;
    template <class T>
    using SequenceType = std::vector<AOG_LIB::Symbolic_State<T> >;
    using VertexId = unsigned int;

    template<class StateType>
    class rule
    {
    public:
        typedef std::shared_ptr<rule> ptr;
        typedef std::shared_ptr<const rule> const_ptr;
        rule(const rule&) = delete;
        rule(){}
        rule& operator=(const rule&) = delete;
        // operator bool()
        // {
        //     return succeeded_;
        // }
        size_t size(unsigned int right) const
        {
            return right_[right].size();
        }
    private:
        Symbolic_State<StateType> left_;
        std::vector<std::vector<Symbolic_State<StateType> > > right_;
        // bool succeeded_;

        //  TODO: check whether "Symbolic_State<StateType>" or vector"Symbolic_State<StateType>"
        rule(const Symbolic_State<StateType>& left, const Symbolic_State<StateType>& right)
            : left_(left)
        {
            right_.emplace_back(1, right);
        }

        // TODO: check constructor
        rule(const Symbolic_Rule<StateType> &str)
        {
            left_ = str.GetSource();

            right_.push_back(str.GetResults());
        }
    public:
        // make right a single state vector
        static ptr create(const Symbolic_State<StateType> & left, const Symbolic_State<StateType>& right)
        {
            return ptr(new rule(left, right));
        }
        static ptr create(const Symbolic_Rule<StateType> & line)
        {
            return ptr(new rule(line));
        }
        const Symbolic_State<StateType>& left() const
        {
            return left_;
        }
        const std::vector<std::vector<Symbolic_State<StateType> > >& right() const
        {
            return right_;
        }
        friend std::ostream& operator <<(std::ostream& os, const rule& r)
        {
            os << r.left_.GetContent() << " -> ";
            unsigned int i = 0;
            for (const std::vector<Symbolic_State<StateType> >& alternative : r.right_) {
                for (const Symbolic_State<StateType> & symbol: alternative) {
                    os << symbol.GetContent() << " ";
                }
                if (i++ < r.right_.size() - 1) {
                    os << "| ";
                }
            }
            return os;
        }
    };
    
    template<class StateType>
    class grammar
    {
    public:
        // Create from a stream
        // grammar(std::istream& is)
        grammar(){}
        grammar(const std::vector<Symbolic_Rule<StateType> >& rules, const std::vector<Symbolic_State<StateType> >& top_level_rules)
        {
            tops_ = top_level_rules;
            for (int i = 0; i < rules.size(); i++){
                typename rule<StateType>::ptr r = rule<StateType>::create(rules[i]);
                rules_.push_back(r);
            }
            // Get the terminals
            // TODO: debug here, changed set to unordered_set
            std::unordered_set<Symbolic_State<StateType> > nonterminals;
            std::unordered_set<Symbolic_State<StateType> > symbols;
            for (typename rule<StateType>::const_ptr r : rules_) {
                nonterminals.insert(r->left());
                for (const std::vector<Symbolic_State<StateType> >& alternative : r->right()) {
                    for (const Symbolic_State<StateType> & symbol : alternative) {
                        symbols.insert(symbol);
                    }
                }
            }
            for (const Symbolic_State<StateType> & symbol : symbols) {
                if (nonterminals.find(symbol) == nonterminals.end()) {
                    terminals_.push_back(symbol);
                }
            }
            // std::cerr << "private rules size inside grammar is: " << this->rules().size() << std::endl;
        }

        // grammar constructor used when doing delete;
        // The sequence could have and node or or node as terminals
        grammar(const std::vector<Symbolic_Rule<StateType>> &rules, const std::vector<Symbolic_State<StateType>> &top_level_rules, 
                const std::vector<Symbolic_State<StateType>> &added_terminals)
        {
            tops_ = top_level_rules;

            for (int i = 0; i < rules.size(); i++)
            {
                typename rule<StateType>::ptr r = rule<StateType>::create(rules[i]);
                // std::cerr << "Grammar initialization: rules.left is: " << r->left().GetContent() << std::endl;
                rules_.push_back(r);
            }
            // Get the terminals
            // TODO: debug here, changed set to unordered_set
            std::unordered_set<Symbolic_State<StateType>> nonterminals;
            std::unordered_set<Symbolic_State<StateType>> symbols;
            for (typename rule<StateType>::const_ptr r : rules_)
            {
                nonterminals.insert(r->left());
                for (const std::vector<Symbolic_State<StateType>> &alternative : r->right())
                {
                    for (const Symbolic_State<StateType> &symbol : alternative)
                    {
                        symbols.insert(symbol);
                    }
                }
            }
            for (const Symbolic_State<StateType> &symbol : symbols)
            {
                if (nonterminals.find(symbol) == nonterminals.end())
                {
                    terminals_.push_back(symbol);
                }
            }

            for (const Symbolic_State<StateType> added_terminal : added_terminals)
            {
                if (std::find(terminals_.begin(), terminals_.end(), added_terminal) == terminals_.end() && !added_terminal.GetIsBasic())
                    terminals_.push_back(added_terminal);
            }

            // std::cerr << "Print terminals at construction of grammar:\n";
            // print_terminal();
        }



        void update(const std::vector<Symbolic_Rule<StateType> >& rules, const std::vector<Symbolic_State<StateType> >& top_level_rules, const SequenceType<StateType> &input_seq)
        {
            this->rules_.clear();
            this->start_rules_.clear();
            this->terminals_.clear();
            this->tops_.clear();
            /** std::cerr << "Update Grammar" << std::endl; */
            tops_ = top_level_rules;

            for (int i = 0; i < rules.size(); i++){
                typename rule<StateType>::ptr r = rule<StateType>::create(rules[i]);
                rules_.push_back(r);
            }
            // Get the terminals
            // TODO: debug here, changed set to unordered_set
            std::unordered_set<Symbolic_State<StateType> > nonterminals;
            std::unordered_set<Symbolic_State<StateType> > symbols;
            for (typename rule<StateType>::const_ptr r : rules_) {
                nonterminals.insert(r->left());
                for (const std::vector<Symbolic_State<StateType> >& alternative : r->right()) {
                    for (const Symbolic_State<StateType> & symbol : alternative) {
                        symbols.insert(symbol);
                    }
                }
            }
            for (const Symbolic_State<StateType> & symbol : symbols) {
                if (nonterminals.find(symbol) == nonterminals.end()) {
                    terminals_.push_back(symbol);
                }
            }

            //add and/or nodes in the input_seq to the terminals_
            for(auto iter = input_seq.begin(); iter != input_seq.end(); iter++)
            {
                if(!(*iter).GetIsBasic())
                {
                    if(std::find(terminals_.begin(), terminals_.end(), (*iter)) == terminals_.end())
                        terminals_.push_back((*iter));
                }
            }
        }

        const std::vector<typename rule<StateType>::ptr>& rules() const
        {
            return rules_;
        }
    
        // Get all of the rules where symbol is the subject
        template <class OutputIt>
        OutputIt get_rules_for_left(const Symbolic_State<StateType> & symbol, OutputIt it) const
        {
            for (typename rule<StateType>::const_ptr r : rules_) {
                if (r->left() == symbol) {
                    *it++ = r;
                }
            }
            return it;
        }
    
        // Is this symbol a terminal (doesn't occur as the subject of a rule)?
        bool symbol_is_terminal(const Symbolic_State<StateType> & symbol) const
        {
            return std::find(terminals_.begin(), terminals_.end(), symbol) != terminals_.end();
        }

        void update_terminal(const Symbolic_State<StateType>& symbol)
        {
            terminals_.push_back(symbol);
        }

            // Get the rule(s) whose left-hand side is the start symbol
            template <class OutputIt>
            OutputIt get_start_rules(OutputIt it) const
        {
            Symbolic_State<StateType> start_symbol;
            bool started = false;
            for (typename rule<StateType>::const_ptr r : rules_) {
                if (!started || r->left() == start_symbol) {
                    *it++ = r;
                }
                if (!started) {
                    started = true;
                    start_symbol = r->left();
                }
            }
            return it;
        }
    
        // Get the rules where this symbol occurs in a right hand side alternative
        template <class OutputIt>
        OutputIt get_rules_for_symbol(const Symbolic_State<StateType> & symbol, OutputIt it) const
        {
            for (typename rule<StateType>::const_ptr r : rules_) {
                for (const std::vector<Symbolic_State<StateType> >& alternative : r->right()) {
                    if (std::find(alternative.begin(), alternative.end(), symbol) != alternative.end()) {
                        *it++ = r;
                    }
                }
            }
            return it;
        }
    
        // Get the subject of the start symbol
        const Symbolic_State<StateType> & get_start_left() const
        {
            return rules_.front()->left();
        }
    
        // Is this symbol the left (subject) of a terminal?
        bool symbol_is_left_of_terminal(const Symbolic_State<StateType> & symbol) const
        {
            for (typename rule<StateType>::const_ptr r : rules_) {
                if (r->left() == symbol) {
                    for (const std::vector<Symbolic_State<StateType> >& alternative : r->right()) {
                        if (alternative.size() == 1 && symbol_is_terminal(alternative[0])) {
                            return true;
                        }
                    }
                }
            }
            return false;
        }
    
        const std::vector<Symbolic_State<StateType> > get_top_level_rules() const
        { return tops_; }

        // Pretty-print
        friend std::ostream& operator <<(std::ostream& os, const grammar& gram)
        {
            for (typename rule<StateType>::const_ptr r : gram.rules_) {
                os << *r;
                os << '\n';
            }
            return os;
        }

        void print_terminal() const
        {
            for(int i = 0; i < terminals_.size(); i++)
            {
                std::cerr << "(" << terminals_[i].GetContent() << ", " << terminals_[i].GetId() << ") ";
            }
        }

      private:
        std::vector<Symbolic_State<StateType>> terminals_;
        std::vector<typename rule<StateType>::ptr> rules_;
        std::vector<typename rule<StateType>::ptr> start_rules_;
        std::vector<Symbolic_State<StateType> > tops_;
    };

    template<class StateType>
    struct state
    {
        typename rule<StateType>::const_ptr rule_; // The grammar rule
        unsigned int right_; // Index of right hand side alternative
        unsigned int dot_; // Position of dot within symbols on right hand side
        unsigned int i_, j_; // Positions within the input
        char added_by_; // Which function added this state
        std::vector<std::vector<std::pair<int, int> > > back_pointer_; // back pointer to next level parsed rules: statelist #, position in statelist
        state(typename rule<StateType>::const_ptr rule, unsigned int right, unsigned int i, unsigned int j)
            : rule_(rule), right_(right), dot_(0), i_(i), j_(j), added_by_(0)
        {        }
    
        // Is dot all the way to the right?
        bool completed() const
        {
            return dot_ == rule_->right()[right_].size();
        }
    
        // Get symbol to the right of dot
        Symbolic_State<StateType> next_symbol() const
        {
            return rule_->right()[right_][dot_];
        }

    };

    // Pretty-print state
    template <class StateType>
    std::ostream& operator <<(std::ostream& os, const state<StateType>& st)
    {
        const std::vector<Symbolic_State<StateType> >& right = st.rule_->right()[st.right_];
        size_t rlen = right.size();
        os << '(';
        os << '('<<st.rule_->left().GetContent()<<", "<<st.rule_->left().GetId()<<")" << " -> ";
        unsigned int s;
        for (s = 0; s < rlen; ++s) {
            if (s == st.dot_) {
                os << "@ ";
            }
            os << "("<<right[s].GetContent()<<", "<<right[s].GetId()<<")";
            if (s < rlen - 1) {
                os << ' ';
            }
        }
        if (s == st.dot_) {
            os << " @";
        }
        os << ", [" << st.i_ << " , " << st.j_ << "]) ";
        switch (st.added_by_) {
            case 'P':
                os << "predictor";
                break;
            case 'S':
                os << "scanner";
                break;
            case 'C':
                os << "completer";
                break;
            default:
                os << "start state";
        }
        return os;
    }
        
    template <class StateType>
    bool operator == (const rule<StateType>& rule1, const rule<StateType>& rule2){
        if(rule1.left() == rule2.left() && rule1.right() == rule2.right()){
            return true;
        }
        else{
            return false;
        }
    }
        
    // Needed to check for duplicate states
    template <class StateType>
    bool operator ==(const state<StateType>& state1, const state<StateType>& state2)
    {
        if (*state1.rule_ == *state2.rule_
            && state1.right_ == state2.right_
            && state1.dot_ == state2.dot_
            && state1.i_ == state2.i_
            && state1.j_ == state2.j_)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    
    // A statelist is a list of states
    template <class StateType>
    struct StateList
    {
        typedef std::vector<state<StateType> > statelist;
    };
    // A chart is a vector of statelists
    template <class StateType>
    struct chart
    {
        const grammar<StateType>& grammar_;
        std::vector<typename StateList<StateType>::statelist> chart_;
        std::vector<typename std::unordered_map<state<StateType>, unsigned> > chart_map_;
        std::vector<unsigned int> succeeded_chart_idx;
    
        chart(const grammar<StateType> & grammar)
            : grammar_(grammar)
        {
            typename StateList<StateType>::statelist top_list;
            std::unordered_map<state<StateType>, unsigned > top_list_map;
            // std::cerr << "Initialize grammar initialization:before get_top_level_rules" << std::endl;            
            std::vector<Symbolic_State<StateType> > top_level_rules = grammar.get_top_level_rules();
            // std::cerr << "Initialize grammar initialization:after get_top_level_rules" << std::endl;                        
            unsigned counter = 0;                   
            for (Symbolic_State<StateType> source : top_level_rules)
            {
                state<StateType> source_state(rule<StateType>::create(Symbolic_State<StateType>("$", false), source), 0, 0, 0);
                top_list.push_back(source_state);
                top_list_map.emplace(source_state, counter);
                counter++;
            }
            // std::cerr << "Initialize grammar initialization:before push_back" << std::endl;
            chart_.push_back(top_list);
            chart_map_.push_back(top_list_map);
        }
    
        // Add state st to statelist s
        void add_state(state<StateType> & st, unsigned int s)
        {
            if (st.rule_->left().GetContent() == "$")
                succeeded_chart_idx.push_back(s);
            if (s < chart_.size()) {
                // Adding to the last statelist
                typename std::unordered_map<state<StateType>, unsigned>::const_iterator it = chart_map_[s].find(st);
                if(it == chart_map_[s].end()){
                    size_t index = chart_[s].size();
                    chart_[s].push_back(st);
                    st.back_pointer_.clear();
                    chart_map_[s].emplace(st, index);
                }
                else if(!st.back_pointer_.empty()){ // rule out states added by predictors
                    state<StateType>& st_stored = chart_[s][it->second];
                    assert(!st_stored.back_pointer_.empty());
                    assert(st.back_pointer_[0].size() == st_stored.back_pointer_[0].size());
                    for(const auto& bp : st.back_pointer_){
                        st_stored.back_pointer_.push_back(bp);
                    }
                }
            }
            else {
                // Adding to a new statelist
                chart_.emplace_back(1, st);
                st.back_pointer_.clear();
                chart_map_.emplace_back(typename std::unordered_map<state<StateType>, unsigned>({{st, 0}}));
            }
        }
    
        // Add predictions for the next symbol in this state
        void predictor(const state<StateType>& st)
        {   
            std::vector<typename rule<StateType>::const_ptr> rules;
            grammar_.get_rules_for_left(st.next_symbol(), std::back_inserter(rules));
            for (typename rule<StateType>::const_ptr r : rules) {
                for (unsigned int a = 0; a < r->right().size(); ++a) {
                    state<StateType> prediction = state<StateType>(r, a, st.j_, st.j_);
                    prediction.added_by_ = 'P';
                    add_state(prediction, st.j_);
                }
            }
        }
    
        // Scan input for next symbol
        void scanner(const state<StateType>& st, const std::vector<Symbolic_State<StateType> >& input, bool back_track = true)
        {
            // position of scanned word
            const Symbolic_State<StateType> & word = input[st.j_];
            if (word == st.rule_->right()[st.right_][st.dot_]) {
                state<StateType> scanned = state<StateType>(st.rule_, st.right_, st.i_, st.j_ + 1);
                scanned.dot_ = st.dot_ + 1;
                scanned.added_by_ = 'S';
                if(back_track && st.back_pointer_.size() > 0){
                    scanned.back_pointer_ = st.back_pointer_;
                }
                add_state(scanned, st.j_ + 1);
            }
        }
    
        // Complete states
        void completer(const state<StateType>& param_st, const int chart_i, const int chart_j, bool back_track = true)
        {
            state<StateType> st = param_st;
            Symbolic_State<StateType> left = st.rule_->left();
            std::vector<Symbolic_State<StateType> > right = st.rule_->right()[st.right_];
            const unsigned int i = st.i_;
            const unsigned int j = st.j_;
            
            unsigned kk = chart_[0][0].rule_->right().size();
            
            int k = 0;
            for (const typename StateList<StateType>::statelist& list = chart_[st.i_]; k < list.size(); k++) {
                if (list[k].j_ == i
                        && !list[k].completed()
                        && list[k].next_symbol() == st.rule_->left())
                {
                    Symbolic_State<StateType> left = list[k].rule_->left();
                    std::vector<Symbolic_State<StateType> > right = list[k].rule_->right()[st.right_];
                    state<StateType> completed = state<StateType>(list[k].rule_, list[k].right_, list[k].i_, j);
                    completed.dot_ = list[k].dot_ + 1;
                    completed.added_by_ = 'C';
                    if (back_track){
                        if(list[k].back_pointer_.size() != 0){
                            completed.back_pointer_ = list[k].back_pointer_; 
                        }
                        if(completed.back_pointer_.size() == 0){
                            completed.back_pointer_.emplace_back(1, std::pair<int, int> (chart_i, chart_j));
                        }
                        else{
                            for(auto& bp : completed.back_pointer_){
                                bp.emplace_back(chart_i, chart_j);
                            }
                        }
                    }
                    add_state(completed, j);
                }
            }
        }
    
        int parsed_position() const
        {
            // chart 0 is dummy
            if(succeeded_chart_idx.empty())
                return 0;
            
            return succeeded_chart_idx.back();
        }

        // The main algorithm
        int parse(const std::vector<Symbolic_State<StateType> >& input, std::ostream *os,
                  std::vector<std::vector<typename StateList<StateType>::statelist> >& parsed_results, bool &parsing_success, bool back_track=true)
        {
            parsing_success = false;

            for (unsigned int i = 0; i <= input.size(); ++i) {
                if (chart_.size() > i) { // Check for running out of statelists when parse fails
                    typename StateList<StateType>::statelist::iterator it = chart_[i].begin();
                    int pos = 0;
                    while(it != chart_[i].end())
                    {
                        const state<StateType> st = *it;
                            
                        if (!st.completed() && !grammar_.symbol_is_terminal(st.next_symbol())) {
                            predictor(st);
                        }
                        else if (!st.completed()) {
                            if (i < input.size()) {
                                scanner(st, input, back_track);
                            }
                        }
                        else {
                            completer(st, i, pos, back_track);
                        }
                        pos++;
                        it = chart_[i].begin()+pos;
                    }
                    // if (os) {
                    //     *os << *this;
                    //     *os << '\n';
                    // }
                }
            }

            int parsed_pos = parsed_position();
            if (parsed_pos != 0 && !succeeded_chart_idx.empty())
            {
                parsing_success = true;
                std::vector<typename StateList<StateType>::statelist> parsed_chart(chart_.begin(),
                                                                    chart_.begin() + parsed_pos + 1);
                parsed_results.push_back(parsed_chart);
                return parsed_pos;
            }
            return 0;
        }
    
        // Pretty-print chart
        friend std::ostream& operator <<(std::ostream& os, const chart &ch)
        {
            for (unsigned int i = 0; i < ch.chart_.size(); ++i) {
                os << "S" << i << ": ";
                os << '[';
                unsigned int j = 0;
                for (const state<StateType>& st : ch.chart_[i]) {
                    os << st;
                    if (j < ch.chart_[i].size() - 1) {
                        os << ",\n";
                    }
                    ++j;
                }
                os << ']';
                os << "\n\n";
            }
            return os;
        }
    };

    template <class StateType, class AttributeType>
    class EarleyParser
    {
    public:
        EarleyParser(){}
        EarleyParser(const grammar<StateType>& grammar)
            : grammar_(grammar)
        {}
        template <class InputIt>
        int parse(InputIt begin, InputIt end)
        {
            std::vector<Symbolic_State<StateType> > input;
            std::copy(begin, end, std::back_inserter(input));
            return chart<StateType>(grammar_).parse(input, nullptr, parsed_results);
        }
        template <class InputIt>
        int parse(InputIt begin, InputIt end, std::ostream& os, bool & parsing_success, bool back_track = true)
        {
            std::vector<Symbolic_State<StateType> > input;
            std::copy(begin, end, std::back_inserter(input));
            // std::cerr << "EarleyParser: before chart construction " << std::endl;
            int pos = chart<StateType>(grammar_).parse(input, &os, parsed_results, parsing_success, back_track);
            // std::cerr << "EarleyParser: position is: " << pos << std::endl;
            return pos;
        }

        double prob(std::shared_ptr<AOG<StateType, AttributeType> > graph_ptr){


            // whole parse tree integrated from sub parse trees
            // not suitable due to ways to parse nonterminals on the right side are related to each other (nonterminals)
            double total_prob = 0;
            assert(this->parsed_results.size() == 1);

            for (const std::vector<typename StateList<StateType>::statelist> &chart : this->parsed_results)
            {
                bool found_root = false;
                int last_statelist_index = chart.size() - 1;
                for (const state<StateType>& st_last : chart[last_statelist_index]){

                    // find top most root
                    if (st_last.i_ == 0 && st_last.j_ != 0 && st_last.rule_->left().GetContent() == "$" && st_last.completed()){
                        //TODO
                        assert(!found_root);
                        found_root = true;

                        std::unordered_map<std::pair<int, int>, double, boost::hash<std::pair<int, int> > > m;

                        std::function<double (const std::pair<int, int>&)> prob_recursive = [&](const std::pair<int,int>& pair) -> double {
                            const state<StateType>& st = chart[pair.first][pair.second];

                            if(m.find(pair) != m.end()){
                                return m[pair];
                            }
                            else{
                                double or_prob = 1;
                                VertexId left = graph_ptr->GetVertexIdByState(st.rule_->left());
                                if (!graph_ptr->GetVertexContent(left)->IsAnd()){
                                    std::unordered_map<VertexId, double> outEdgeWeights = graph_ptr->GetOutEdgeWeights(left, true);
                                    std::vector<VertexId> right_hand_state_ids;
                                    for (const auto& right_hand_state : st.rule_->right()[st.right_]){
                                        right_hand_state_ids.push_back(graph_ptr->GetVertexIdByState(right_hand_state));
                                    }
                                    std::vector<VertexId> dummy_vertices = graph_ptr->ChildrenVertices(left);

                                    bool found_destination = false;
                                    for (VertexId id : dummy_vertices)
                                    {
                                        std::vector<VertexId> children_vertices = graph_ptr->ChildrenVertices(id);
                                        
                                        if (children_vertices == right_hand_state_ids)
                                        {
                                            found_destination = true;
                                            or_prob = outEdgeWeights[id];
                                            break;
                                        }
                                    }
                                    assert(found_destination);
                                }

                                std::unordered_set<std::vector<std::pair<int, int> >, container_hash<std::vector<std::pair<int, int> > > > back_pointer_set;

                                double children_prob = 0;
                                for (const auto& bp_vector : st.back_pointer_){
                                    back_pointer_set.insert(bp_vector);
                                    double children_prob_inner = 1;
                                    for(const auto& bp : bp_vector){
                                        children_prob_inner *= prob_recursive(bp);
                                    }
                                    children_prob += children_prob_inner;
                                }

                                if(st.back_pointer_.empty()){
                                    for(const Symbolic_State<StateType>& state : st.rule_->right()[st.right_]){
                                        assert(state.GetIsBasic());
                                    }
                                    children_prob = 1;
                                }
                                assert(back_pointer_set.size() == st.back_pointer_.size());

                                double pair_prob = or_prob * children_prob;
                                m.emplace(pair, pair_prob);

                                return pair_prob;
                            }

                        };

                        for(const auto& bp_vector : st_last.back_pointer_){
                            assert(bp_vector.size() == 1);

                            total_prob += prob_recursive(bp_vector[0]);

                            assert(!m.empty());
                        }

                        break;
                    }
                }
            }

            return total_prob;
        }

        std::vector<std::vector<typename StateList<StateType>::statelist > > GetPartialParse(){
            return parsed_results;
        }
        grammar<StateType> GetGrammar(){
            return grammar_;
        }

        void update(const std::vector<Symbolic_Rule<StateType> > & rules, const SequenceType<StateType> &input_seq)
        {
            SequenceType<StateType> top_level_rules = this->get_top_level_rules(rules);
            this->grammar_.update(rules, top_level_rules, input_seq);
        }

        std::vector<Symbolic_State<StateType> > get_top_level_rules(std::vector<Symbolic_Rule<StateType> > rules)
        {
            std::vector<Symbolic_State<StateType> > top_level_rules;
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
                    top_level_rules.push_back(source);
            }
            
            return top_level_rules;
        }
    private:
        grammar<StateType> grammar_;
        std::vector<std::vector<typename StateList<StateType>::statelist > > parsed_results;
    };

}

namespace std
{

    template <typename StateType>
    size_t hash_value(const AOG_LIB_UTIL::rule<StateType> &rule){
        std::size_t seed = 0;
        boost::hash_combine(seed, rule.left());
        boost::hash_combine(seed, rule.right());
        return seed;
    }

    template <typename StateType>
    struct hash<AOG_LIB_UTIL::rule<StateType> >
    {
        size_t operator()(const AOG_LIB_UTIL::rule<StateType> &rule) const noexcept
        {
            boost::hash<AOG_LIB_UTIL::rule<StateType> > hasher;
            return hasher(rule);
        }
    };

    template <typename StateType>
    struct hash<AOG_LIB_UTIL::state<StateType> >
    {
        size_t operator()(const AOG_LIB_UTIL::state<StateType> &state) const noexcept
        {
            std::size_t seed = 0;
            boost::hash_combine(seed, *state.rule_);
            boost::hash_combine(seed, state.right_);
            boost::hash_combine(seed, state.dot_);
            boost::hash_combine(seed, state.i_);
            boost::hash_combine(seed, state.j_);
            // boost::hash_combine(seed, state.back_pointer_);
            return seed;
        }
    };
}


#endif //AOG_LIB_EARLEY_EVALUATION_H