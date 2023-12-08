# model.py
# Python version of WatDiv
# T. Masuda, 2023/12/1

import random
import sys
import math
import numpy as np
from enum import Enum, auto
import bisect
from datetime import datetime

# #include "dictionary.h"
from src.dictionary import Dictionary
# #include "model.h"
# #include "statistics.h"
# from src.statistics import Statistics  # statistics is merged into this file.
# #include "volatility_gen.h"
from src.volatility_gen import VolatilityGen

# #include <chrono>
# #include <fstream>

# #include <iostream>

# #include <random>
# #include <set>
# #include <map>
# #include <unordered_set>
# #include <sstream>
# #include <stack>

# #include <math.h>

# #include <boost/algorithm/string.hpp>
# #include <boost/algorithm/string/predicate.hpp>
# #include <boost/algorithm/string/split.hpp>
# #include <boost/date_time/local_time/local_time.hpp>
# #include "boost/date_time/posix_time/posix_time.hpp"
# #include "boost/date_time/posix_time/posix_time_types.hpp"
# #include <boost/lexical_cast.hpp>

# namespace bpt = boost::posix_time;

# -------------------------------------------------------------------------------------------------------
# statistics.py
# Python version of WatDiv/statistics
# T. Masuda, 2023/12/1

# import random
# from enum import Enum, auto
# #include "statistics.h"

# #include <algorithm>
# #include <stdlib.h>     /* srand, rand */
# #include <time.h>       /* time */

# #include <boost/algorithm/string/replace.hpp>
# #include <boost/lexical_cast.hpp>

# namespace QUERY_STRUCTURE {
#     enum enum_t {
#         PATH,
#         STAR,
#         SNOWFLAKE,
#         COMPLEX,
#         UNDEFINED
#     };
# };


class QueryStructure(Enum):
    PATH = auto
    STAR = auto
    SNOWFLAKE = auto
    COMPLEX = auto
    UNDEFINED = auto

# namespace QUERY_CATEGORY {
#     enum enum_t {
#         LOW_SELECTIVITY,
#         MEDIUM_SELECTIVITY,
#         HIGH_SELECTIVITY,
#         UNDEFINED
#     };
# };


class QueryCategory(Enum):
    LOW_SELECTIVITY = auto
    MEDIUM_SELECTIVITY = auto
    HIGH_SELECTIVITY = auto
    UNDEFINED = auto

# ostream& operator<<(ostream& os, const QUERY_STRUCTURE::enum_t & query_structure){
#     switch (query_structure){
#         case QUERY_STRUCTURE::PATH:{
#             os << "[path]";
#             break;
#         }
#         case QUERY_STRUCTURE::STAR:{
#             os << "[star]";
#             break;
#         }
#         case QUERY_STRUCTURE::SNOWFLAKE:{
#             os << "[snowflake]";
#             break;
#         }
#         case QUERY_STRUCTURE::COMPLEX:{
#             os << "[complex]";
#             break;
#         }
#         case QUERY_STRUCTURE::UNDEFINED:{
#             os << "[undefined]";
#             break;
#         }
#     }
#     return os;
# }

# ostream& operator<<(ostream& os, const QUERY_CATEGORY::enum_t & query_category){
#     switch (query_category){
#         case QUERY_CATEGORY::LOW_SELECTIVITY:{
#             os << "[low]";
#             break;
#         }
#         case QUERY_CATEGORY::MEDIUM_SELECTIVITY:{
#             os << "[medium]";
#             break;
#         }
#         case QUERY_CATEGORY::HIGH_SELECTIVITY:{
#             os << "[high]";
#             break;
#         }
#         case QUERY_CATEGORY::UNDEFINED:{
#             os << "[undefined]";
#             break;
#         }
#     }
#     return os;
# }

# ostream& operator<<(ostream& os, const statistics_st & stats){
#     os <<"\t" << "(" << stats._vertex1 << ", " << stats._edge << ", " << stats._vertex2 << ", " << stats._right_distribution << ")" << "\n";
#     return os;
# }


class StatisticsST:
    """

    Attributes:
        _vertex1 (str): subject of a triple
        _vertex2 (str): object of a triple
        _edge (str): predicate of a triple
        _right_distribution (DistributionTypes):
        _group_id (int):
        _pr_left (float)
        _pr_right (float)
    """
#     string _vertex1;
#     string _vertex2;
#     string _edge;
#     DISTRIBUTION_TYPES::enum_t _right_distribution;
#     int _group_id;
#     double _pr_left;
#     double _pr_right;


# statistics_st::statistics_st (const string & vertex1, const string & vertex2, const string & edge, DISTRIBUTION_TYPES::enum_t right_distribution, int group_id){
    def __init__(self, vertex1: str, vertex2: str, edge: str, right_distribution, group_id: int):
        """
        Constructor of StatisticsST

        Args:
            vertex1 (str):

        Returns:

        """
#     _vertex1 = vertex1;
        self._vertex1 = vertex1
#     _vertex2 = vertex2;
        self._vertex2 = vertex2
#     _edge = edge;
        self._edge = edge
#     _right_distribution = right_distribution;
        self._right_distribution = right_distribution
#     _group_id = group_id;
        self._group_id = group_id
#     _pr_left = 0.0;
        self._pr_left = 0.0
#     _pr_right = 0.0;
        self._pr_right = 0.0
# }
        pass  # end of __init__

# bool statistics_st::operator< (const statistics_st & rhs) const{
#     if (_edge.compare(rhs._edge)!=0){
#         return _edge.compare(rhs._edge)<0;
#     }
#     if (_vertex1.compare(rhs._vertex1)!=0){
#         return _vertex1.compare(rhs._vertex1)<0;
#     }
#     if (_vertex2.compare(rhs._vertex2)!=0){
#         return _vertex2.compare(rhs._vertex2)<0;
#     }
#     return false;
# }
    pass  # end of StatisticsST


class Statistics:
    """

    Attributes:
        graph (dict[str, set[StatisticsST]]):
        _spo_index (list[TripleST]):
        _ops_index  (list[TripleST]):
        _statistics_table (dic[str, (float, float)]
        _model (Model)
    """

# statistics::statistics(const model * mdl, const vector<triple_st> & triple_array, int maxQSize, int qCount, int constCount, bool constJoinVertexAllowed, bool dupEdgesAllowed){
    def __init__(self, mdl, triple_array, maxQSize, qCount, constCount, constJoinVertexAllowed, dupEdgesAllowed):
        """

        Args:

        Returns:

        """
# const model * _model;
#         map<string, set<statistics_st> > graph; // @suppress("Invalid template argument")
        self.graph = {}
#         vector<triple_st> _spo_index;
        self._spo_index = []
#         vector<triple_st> _ops_index;
        self._ops_index = []
#         map<string, pair<double, double> > _statistics_table; // @suppress("Invalid template argument")
        self._statistics_table = {}
#     srand (time(NULL));
#     _model = mdl;
        self._model = mdl
#     index_triples(triple_array);
        self.index_triples(triple_array)
#     extract_schema(*_model);
        self.extract_schema(self._model)
#     //print_graph();

#     ///compute();
#     //print_stats();

#     map<pair<QUERY_STRUCTURE::enum_t, QUERY_CATEGORY::enum_t>, map<string, string> > query_map;
        query_map = {}

        i = 0
#     for (int i=0; i<qCount; ){
        while i < qCount:  # repeat until the number of success reaches qCount
#         int qSize = (rand() % maxQSize) + 1;
            qSize = int(random.uniform(0, maxQSize-1)) + 1  # number of triples in a query
#         if (traverse_graph(qSize, constCount, constJoinVertexAllowed, dupEdgesAllowed, query_map)){
            if self.traverse_graph(qSize, constCount, constJoinVertexAllowed, dupEdgesAllowed, query_map):
#             i++;
                i += 1  # succeeded
#         }  // end of if
#     }  // end of for

# #     /*
# #     // Randomly select a subset of queries from each category...
# #     cout << "Randomly selected query templates ..." << "\n";
#         print("Randomly selected query templates ...")
# #     map<pair<QUERY_STRUCTURE::enum_t, QUERY_CATEGORY::enum_t>, map<string, string> >::iterator itr1 = query_map.begin();
# #     for (; itr1!=query_map.end(); itr1++){
#         for itr1_key, itr1_value in query_map.items():
# #         if (itr1->first.first != QUERY_STRUCTURE::UNDEFINED && itr1->first.second != QUERY_CATEGORY::UNDEFINED){
#             if itr1_key[0] != QueryStructure.UNDEFINED and itr1_key[1] != QueryCategory.UNDEFINED:
# #             cout << "CATEGORY ::: " << itr1->first.first << " - " << itr1->first.second << "\n";
#                 print(f'CATEGORY ::: {itr1_key[0]} - {itr1_key[1]}')
# #             map<string, string> query_description_map = itr1->second;
#                 query_description_map = itr1_value
# #             vector<string> keys;
#                 keys = []
# #             for (map<string, string>::iterator itr2=query_description_map.begin(); itr2!=query_description_map.end(); itr2++){
#                 for itr2 in query_description_map:
# #                 keys.push_back(itr2->first);
#                     keys.append(itr2[0])
# #             }
# #             for (int k=0; k<5; k++){
#                 for k in range(5):
#                     xxx = random.uniform(0, len(keys) - 1)
# #                 string random_query = keys[rand()%keys.size()];
#                     random_query = keys[xxx]
# #                 string description = query_description_map[random_query];
#                     description = query_description_map[random_query]
# #                 cout << random_query << description << "\n";
#                     print(f'{random_query} {description}')
# #             }
# #         }
# #     }
#             pass  # end of for
# #     */
# }
        pass  # end of statistics::statistics(__init__)

# statistics::~statistics(){
# }

# void statistics::extract_schema(const model & mdl){
    def extract_schema(self, mdl):
        """

        Args:
            mdl:

        Returns:

        """
#     vector<statistics_st> * result = new vector<statistics_st> ();
        result = []
#     int group_id = 0;
        group_id = 0
#     for (vector<resource_m_t*>::const_iterator itr1=mdl._resource_array.begin(); itr1!=mdl._resource_array.end(); itr1++){
        for itr1 in mdl._resource_array:
#         const resource_m_t * rsrc = *itr1;
            rsrc = itr1
#         for (vector<predicate_group_m_t*>::const_iterator itr2=rsrc->_predicate_group_array.begin(); itr2!=rsrc->_predicate_group_array.end(); itr2++){
            for itr2 in rsrc._predicate_group_array:
#             const predicate_group_m_t * pgroup = *itr2;
                pgroup = itr2
#             for (vector<predicate_m_t*>::const_iterator itr3=pgroup->_predicate_array.begin(); itr3!=pgroup->_predicate_array.end(); itr3++){
                for itr3 in pgroup._predicate_array:
#                 const predicate_m_t * pred = *itr3;
                    pred = itr3
#                 string vertex1 = "", vertex2 = "", edge = "";
                    vertex1 = ''
                    vertex2 = ''
                    edge = ''
#                 vertex1.append(rsrc->_type_prefix);
                    vertex1 += rsrc._type_prefix
#                 if (pgroup->_type_restriction!=NULL){
                    if pgroup._type_restriction is not None:
#                     vertex1.append("@");
                        vertex1 += '@'
#                     vertex1.append(*pgroup->_type_restriction);
                        vertex1 += pgroup._type_restriction
#                 }
#                 switch (pred->_literal_type){
#                     case LITERAL_TYPES::DATE:{
                    if pred._literal_type == LiteralTypes.DATE:
#                         vertex2.append("date");
                        vertex2 += 'data'
#                         break;
#                     }
#                     case LITERAL_TYPES::INTEGER:{
                    elif pred._literal_type == LiteralTypes.INTEGER:
#                         vertex2.append("integer");
                        vertex2 += 'integer'
#                         break;
#                     }
#                     case LITERAL_TYPES::NAME:{
                    elif pred._literal_type == LiteralTypes.NAME:
#                         vertex2.append("name");
                        vertex2 += 'name'
#                         break;
#                     }
#                     case LITERAL_TYPES::STRING:{
                    elif pred._literal_type == LiteralTypes.STRING:
#                         vertex2.append("string");
                        vertex2 += 'string'
#                         break;
#                     }
#                 }
#                 edge.append(pred->_label);
                    edge += pred._label
#                 result->push_back(statistics_st(vertex1, vertex2, edge, pred->_distribution_type, group_id));
                    result.append(StatisticsST(vertex1, vertex2, edge, pred._distribution_type, group_id))
#             }
#             group_id++;
                group_id += 1
#         }
#     }
#     for (vector<association_m_t*>::const_iterator itr1=mdl._association_array.begin(); itr1!=mdl._association_array.end(); itr1++){
        for itr1 in mdl._association_array:
#         const association_m_t * assoc = *itr1;
            assoc = itr1
#         string vertex1 = "", vertex2 = "", edge = "";
            vertex1 = ''
            vertex2 = ''
            edge = ''
#         vertex1.append(assoc->_subject_type);
            vertex1 += assoc._subject_type
#         if (assoc->_subject_type_restriction!=NULL){
            if assoc._subject_type_restriction is not None:
#             vertex1.append("@");
                vertex1 += '@'
#             vertex1.append(*assoc->_subject_type_restriction);
                vertex1 += assoc._subject_type_restriction
#         }
#         vertex2.append(assoc->_object_type);
            vertex2 += assoc._object_type
#         if (assoc->_object_type_restriction!=NULL){
            if assoc._object_type_restriction is not None:
#             vertex2.append("@");
                vertex2 += '@'
#             vertex2.append(*assoc->_object_type_restriction);
                vertex2 += assoc._object_type_restriction
#         }
#         edge.append(assoc->_predicate);
            edge += assoc._predicate
#         result->push_back(statistics_st(vertex1, vertex2, edge, assoc->_right_distribution, group_id));
            result.append(StatisticsST(vertex1, vertex2, edge, assoc._right_distribution, group_id))
#         group_id++;
            group_id += 1
#     }
#     populate_graph(*result);
        self.populate_graph(result)
#     infer_edges();
        self.infer_edges()
#     delete result;
# }
        pass  # end of extract_schema

# void statistics::populate_graph(const vector<statistics_st> & tuples){
    def populate_graph(self, tuples):
        """

        Args:
            tuples:

        Returns:

        """
#     for (vector<statistics_st>::const_iterator itr1=tuples.begin(); itr1!=tuples.end(); itr1++){
        for itr1 in tuples:
#         statistics_st tuple = *itr1;
            tuple = itr1
#         if (graph.find(tuple._vertex1)==graph.end()){
            try:
                xxx = self.graph[tuple._vertex1]
            except KeyError:
#             graph.insert(pair<string, set<statistics_st> >(tuple._vertex1, set<statistics_st>()));
                self.graph[tuple._vertex1] = set()
#         }
#         graph[tuple._vertex1].insert(tuple);
            self.graph[tuple._vertex1].add(tuple)

# //        int pos = string::npos;
# //        if ( (pos = tuple._vertex1.find("@"))!=string::npos){
# //            string inferred_vertex = tuple._vertex1.substr(0, pos);
# //            if (graph.find(inferred_vertex)==graph.end()){
# //                graph.insert(pair<string, set<statistics_st> >(inferred_vertex, set<statistics_st>()));
# //            }
# //            graph[inferred_vertex].insert(tuple);
# //        }

#         if (tuple._vertex2.compare("date") !=0 && tuple._vertex2.compare("integer") !=0 && tuple._vertex2.compare("name") !=0 && tuple._vertex2.compare("string")){
            if tuple._vertex2 != 'date' and tuple._vertex2 != 'integer' and tuple._vertex2 != 'name' and tuple._vertex2 != 'string':
#             if (graph.find(tuple._vertex2)==graph.end()){
                try:
                    xxx = self.graph[tuple._vertex2]
                except KeyError:
#                 graph.insert(pair<string, set<statistics_st> >(tuple._vertex2, set<statistics_st>()));
                    self.graph[tuple._vertex2] = set()
#             }
#             graph[tuple._vertex2].insert(tuple);
                self.graph[tuple._vertex2].add(tuple)

# //            pos = string::npos;
# //            if ( (pos = tuple._vertex2.find("@"))!=string::npos){
# //                string inferred_vertex = tuple._vertex2.substr(0, pos);
# //                if (graph.find(inferred_vertex)==graph.end()){
# //                    graph.insert(pair<string, set<statistics_st> >(inferred_vertex, set<statistics_st>()));
# //                }
# //                graph[inferred_vertex].insert(tuple);
# //            }
#         }
#     }
# }
        pass  # end of

# void statistics::infer_edges(){
    def infer_edges(self):
        """

        Args:

        Returns:

        """
#     for (map<string, set<statistics_st> >::const_iterator itr1=graph.begin(); itr1!=graph.end(); itr1++){
        for itr1_key, iter1_value in self.graph.items():
#         string vertex = itr1->first;
            vertex = itr1_key
#         int pos = string::npos;
            pos = vertex.find('@')
#         if ( (pos = vertex.find("@"))!=string::npos){
            if pos >= 0:
#             string base_vertex = vertex.substr(0, pos);
                base_vertex = vertex[0: pos]
#             // Insert all edges from base class...
#             set<statistics_st> base_edges = graph[base_vertex];
                base_edges = self.graph[base_vertex]
#             for (set<statistics_st>::iterator itr2=base_edges.begin(); itr2!=base_edges.end(); itr2++){
                for itr2 in base_edges:
#                 graph[vertex].insert(*itr2);
                    self.graph[vertex].add(itr2)
#             }
#         }
#     }
# }
        pass  # end of

# bool statistics::traverse_graph(int max_size, int const_count, bool constJoinVertexAllowed, bool dupEdgesAllowed, map<pair<QUERY_STRUCTURE::enum_t, QUERY_CATEGORY::enum_t>, map<string, string> > & query_map) const{
    def traverse_graph(self, max_size, const_count, constJoinVertexAllowed, dupEdgesAllowed, query_map):
        """

        Args:

        Returns:

        """
#     vector<string> v_array;
        v_array = []
#     for (map<string, set<statistics_st> >::const_iterator itr1=graph.begin(); itr1!=graph.end(); itr1++){
        for itr1_key, itr1_value in self.graph.items():
#         v_array.push_back(itr1->first);
            v_array.append(itr1_key)
#     }
        pass  # end of for

#     string v = v_array[rand() % v_array.size()];
        v = v_array[int(random.uniform(0, len(v_array)-1))]

#     set<statistics_st> traversed_edges;
        traversed_edges = set()
#     //multiset<statistics_st> traversed_edges;
#     int itr_counter = 0;
        itr_counter = 0
#     do {
        while True:
#         vector<statistics_st> potential_edges;
            potential_edges = []
#         map<string, set<statistics_st> >::const_iterator fit = graph.find(v);
            fit = self.graph[v]
#         set<statistics_st> incident_edges = fit->second;
            incident_edges = fit
#         for (set<statistics_st>::const_iterator itr1=incident_edges.begin(); itr1!=incident_edges.end(); itr1++){
            for itr1 in incident_edges:
                found = False
                if itr1 in traversed_edges:
                    found = True
#             if (dupEdgesAllowed || traversed_edges.find(*itr1)==traversed_edges.end()){
                if dupEdgesAllowed or not found:
#                 potential_edges.push_back(*itr1);
                    potential_edges.append(itr1)
#             }
#         }
#         if (potential_edges.empty()){
            if len(potential_edges) == 0:
#             break;
                break
#         }

#         statistics_st next_edge = potential_edges[rand() % potential_edges.size()];
            next_edge = potential_edges[int(random.uniform(0, len(potential_edges)-1))]
#         //cout << "\t" << "(" << next_edge._vertex1 << ", " << next_edge._edge << ", " << next_edge._vertex2 << ")" << "\n";
#         traversed_edges.insert(next_edge);
            traversed_edges.add(next_edge)

#         vector<string> potential_next_vertices;
            potential_next_vertices = []
#         if (graph.find(next_edge._vertex1)!=graph.end()){
            try:
                xxx = self.graph[next_edge._vertex1]
#             potential_next_vertices.push_back(next_edge._vertex1);
                potential_next_vertices.append(next_edge._vertex1)
            except KeyError:
                pass
#         }
#         if (graph.find(next_edge._vertex2)!=graph.end()){
            try:
                xxx = self.graph[next_edge._vertex2]
#             potential_next_vertices.push_back(next_edge._vertex2);
                potential_next_vertices.append(next_edge._vertex2)
            except KeyError:
                pass
#         }
            if len(potential_next_vertices) == 0:
                print('potential_next_vertices is empty. ')
                pass  # error
#         v = potential_next_vertices[rand() % potential_next_vertices.size()];
            v = potential_next_vertices[int(random.uniform(0, len(potential_next_vertices)-1))]

#         itr_counter++;
            itr_counter += 1
#     } while (traversed_edges.size()<max_size && itr_counter<(max_size*20));
            if len(traversed_edges) >= max_size:
                break
            if itr_counter >= max_size * 20:
                break
        pass  # end of while

#     // Compute vertices in query graph...
#     set<string> variable_set;
        variable_set = set()
#     string query_str = "";
        query_str = ''
#     int var_count = 0;
        var_count = 0
#     map<string, int> q_vertex_map;
        q_vertex_map = {}
#     map<string, set<string> > variable_map;
        variable_map = {}
#     for (set<statistics_st>::iterator itr1=traversed_edges.begin(); itr1!=traversed_edges.end(); itr1++){
        for itr1 in traversed_edges:
#         string var1 = "", var2 = "";
            var1 = ''
            var2 = ''
#         string v1_base = itr1->_vertex1, v2_base = itr1->_vertex2;
            v1_base = itr1._vertex1
            v2_base = itr1._vertex2
#         int pos = string::npos;
            pos = v1_base.find('@')
#         if ((pos=v1_base.find("@"))!=string::npos){
            if pos >= 0:
#             v1_base = v1_base.substr(0, pos);
                v1_base = v1_base[0: pos]
#         }
            pos = v2_base.find('@')
#         if ((pos=v2_base.find("@"))!=string::npos){
            if pos >= 0:
#             v2_base = v2_base.substr(0, pos);
                v2_base = v2_base[0: pos]
#         }
#         if (q_vertex_map.find(v1_base)==q_vertex_map.end()){
            try:
                xxx = q_vertex_map[v1_base]
            except KeyError:
#             q_vertex_map.insert(pair<string, int>(v1_base, var_count));
                q_vertex_map[v1_base] = var_count
#             var_count++;
                var_count += 1
#         }
            pass  # end of try

#         var1.append("?v");
            var1 += '?v'
#         var1.append(boost::lexical_cast<string>(q_vertex_map[v1_base]));
            var1 += str(q_vertex_map[v1_base])
#         query_str.append("\t");
            query_str += '\t'
#         query_str.append(var1);
            query_str += var1
#         query_str.append("\t");
            query_str += '\t'
#         query_str.append(itr1->_edge);
            query_str += itr1._edge
#         query_str.append("\t");
            query_str += '\t'

#         if (v2_base.compare("date") !=0 && v2_base.compare("integer") !=0 && v2_base.compare("name") !=0 && v2_base.compare("string")){
            if v2_base != 'date' and v2_base != 'integer' and v2_base != 'name' and v2_base != 'string':
#             if (q_vertex_map.find(v2_base)==q_vertex_map.end()){
                try:
                    xxx = q_vertex_map[v2_base]
                except KeyError:
#                 q_vertex_map.insert(pair<string, int>(v2_base, var_count));
                    q_vertex_map[v2_base] = var_count
#                 var_count++;
                    var_count += 1
#             }
                pass  # end of if/try

#             var2.append("?v");
                var2 += '?v'
#             var2.append(boost::lexical_cast<string>(q_vertex_map[v2_base]));
                var2 += str(q_vertex_map[v2_base])
#             query_str.append(var2);
                query_str += var2
#         } else {
            else:
#             var2.append("?v");
                var2 += '?v'
#             var2.append(boost::lexical_cast<string>(var_count));
                var2 += str(var_count)
#             query_str.append(var2);
                query_str += var2
#             var_count++;
                var_count += 1
#         }
            pass  # end of if

#         query_str.append(" ");
            query_str += ' '
#         query_str.append(".");
            query_str += '.'
#         query_str.append("\n");
            query_str += '\n'
#         if (variable_map.find(var1)==variable_map.end()){
            try:
                xxx = variable_map[var1]
            except KeyError:
#             variable_map.insert(pair<string, set<string> >(var1, set<string>()));
                variable_map[var1] = []
#         }
            pass  # end of try
#         variable_map[var1].insert(itr1->_vertex1);
            variable_map[var1].append(itr1._vertex1)
#         if (variable_map.find(var2)==variable_map.end()){
            try:
                xxx = variable_map[var2]
            except KeyError:
#             variable_map.insert(pair<string, set<string> >(var2, set<string>()));
                variable_map[var2] = []
#         }
            pass  # end of try
#         variable_map[var2].insert(itr1->_vertex2);
            variable_map[var2].append(itr1._vertex2)
#         variable_set.insert(var1);
            variable_set.add(var1)
#         variable_set.insert(var2);
            variable_set.add(var2)
#     }
        pass  # end of for

#     /////////////////////////////////////////////////////////////////////////////
#     // Put all traversed edges in a graph...
#     // When constructing the graph, be careful to ignore literals...
#     // For each join vertex (i.e., vertex with more than one incident edge),
#     //  compute join selectivity.
#     /////////////////////////////////////////////////////////////////////////////

#     var_count = 0;
        var_count = 0
#     map<string, set<statistics_st> > query_graph;
        query_graph = {}
#     for (set<statistics_st>::iterator itr1 = traversed_edges.begin(); itr1!= traversed_edges.end(); itr1++){
        for itr1 in traversed_edges:
#         string v1_base = itr1->_vertex1, v2_base = itr1->_vertex2;
            v1_base = itr1._vertex1
            v2_base = itr1._vertex2
#         int pos = string::npos;
            pos = v1_base.find('@')
#         if ((pos=v1_base.find("@"))!=string::npos){
            if pos >= 0:
#             v1_base = v1_base.substr(0, pos);
                v1_base = v1_base[0: pos]
#         }
            pos = v2_base.find('@')
#         if ((pos=v2_base.find("@"))!=string::npos){
            if pos >= 0:
#             v2_base = v2_base.substr(0, pos);
                v2_base = v2_base[0: pos]
#         }
#         if (query_graph.find(v1_base)==query_graph.end()){
            try:
                xxx = query_graph[v1_base]
            except KeyError:
#             query_graph.insert(pair<string, set<statistics_st> >(v1_base, set<statistics_st>()));
                query_graph[v1_base] = set()
#         }
#         query_graph[v1_base].insert(*itr1);
            query_graph[v1_base] = itr1

#         if (v2_base.compare("date") !=0 && v2_base.compare("integer") !=0 && v2_base.compare("name") !=0 && v2_base.compare("string")){
            if v2_base != 'date' and v2_base != 'integer' and v2_base != 'name' and v2_base != 'string':
#             if (query_graph.find(v2_base)==query_graph.end()){
                try:
                    xxx = query_graph[v2_base]
                except KeyError:
#                 query_graph.insert(pair<string, set<statistics_st> >(v2_base, set<statistics_st>()));
                    query_graph[v2_base] = set()
#             }
#             query_graph[v2_base].insert(*itr1);
                    query_graph[v2_base] = itr1
#         } else {
            else:
#             string tmp_id = "?v";
                    tmp_id = '?v'
#             tmp_id.append(boost::lexical_cast<string>(var_count));
                    tmp_id += str(var_count)
#             query_graph[tmp_id].insert(*itr1);
                    query_graph[tmp_id].add(itr1)
#             var_count++;
                    var_count += 1
#         }
            pass  # end of if
#     }
        pass  # end of for

#     /*
#     QUERY_STRUCTURE::enum_t query_structure = QUERY_STRUCTURE::UNDEFINED;
#     QUERY_CATEGORY::enum_t query_category = QUERY_CATEGORY::UNDEFINED;
#     double min_join_selectivity = 1.0;
#     double max_join_selectivity = 1.0;
#     */

#     /*
        """
        #     vector<int> degree_array;
                    degree_array = []
        #     for (map<string, set<statistics_st> >::iterator itr1=query_graph.begin(); itr1!=query_graph.end(); itr1++){
                    for itr1_key, itr1_value in query_graph.items():
        #         const set<statistics_st> edge_list = itr1->second;
                        efge_list = itr1_value
        #         degree_array.push_back(edge_list.size());
                        degree_array.append(len(efge_list))
        #         if (edge_list.size()>=2){
                        if len(efge_list) >= 2:
        #             set<int> correlation_set;
                            correlation_set = set()
        #             for (set<statistics_st>::const_iterator itr2=edge_list.begin(); itr2!=edge_list.end(); itr2++){
                            for itr2 in efge_list:
        #                 if (correlation_set.find(itr2->_group_id)==correlation_set.end()){
                                if itr2._group_id in correlation_set:
                                    pass
                                else:
        #                     correlation_set.insert(itr2->_group_id);
                                    correlation_set.add(itr2._group_id)
        #                     string entity = "";
                                    entity = ''
        #                     bool direction = true;
                                    direction = True
        #                     if (itr2->_vertex1.find(itr1->first)!=string::npos){
                                    if itr1_key in itr2._vertex1:
        #                         entity = itr2->_vertex1;
                                        entity = itr2._vertex1
        #                         direction = true;
                                        direction = True
        #                     } else if (itr2->_vertex2.find(itr1->first)!=string::npos){
                                    elif itr1_key in itr2._vertex2:                            
        #                         entity = itr2->_vertex2;
                                        entity = itr2._vertex2
        #                         direction = false;
                                        direction = False
        #                     } else {
                                    else:
        #                         cout<<"Direction is wrong..."<<"\n";
                                        print("Direction is wrong...")
        #                         cout<<"Vertex "<<itr1->first<<"\n";
                                        print(f'Vertex {itr1_key}')
        #                         cout<<"Edge-1 "<<itr2->_vertex1<<"\n";
                                        print(f'Edge-1 {itr2._vertex1}')
        #                         cout<<"Edge-2 "<<itr2->_vertex2<<"\n";
                                        print(f'Edge-2 {itr2._vertex2}')
        #                         exit(0);
                                        sys.exit(0)
        #                     }
        #                     double max_stats = 0.0;
                                    max_stats = 0.0
        #                     double min_stats = 1.0;
                                    min_stats = 1.0
        #                     for (int dist=0; dist<DISTRIBUTION_TYPES::UNDEFINED; dist++){
                                    for dist in range(DistributionTypes.UNDEFINED):
        #                         string key = get_key(entity, itr2->_edge, direction, (DISTRIBUTION_TYPES::enum_t) dist);
                                        key = self.get_key(entity=entity, predicate=itr2._edge, direction=direction, distribution=dist)
        #                         map<string, pair<double, double> >::const_iterator f_itr = _statistics_table.find(key);
                                        f_itr = self._statistics_table[key]
        #                         pair<double, double> stats = f_itr->second;
                                        stats = f_itr
        #                         //cout<<"\t"<<key<<"\t"<<stats.first<<"\n";
        #                         if (stats.first>max_stats){
                                        if stats[0] > max_stats:
        #                             max_stats = stats.first;
                                            max_stats = stats[0]
        #                         }
        #                         if (stats.first<min_stats){
                                        if stats[0] <min_stats:
        #                             min_stats = stats.first;
                                            min_stats = stats[0]
        #                         }
        #                     }
        #                     min_join_selectivity = min_join_selectivity * min_stats;
                                    min_join_selectivity *= min_stats
        #                     max_join_selectivity = max_join_selectivity * max_stats;
                                    max_join_selectivity *= min_stats
        #                 }
        #             }
        #         }
        #     }
        
        #     sort(degree_array.begin(), degree_array.end());
                        degree_array.sort()
        
        #     if (degree_array.size() >=3 && degree_array[0]==1 && degree_array[1]==1 && degree_array[2]==2 && degree_array[degree_array.size()-1]==2){ // Check if it is a path query...
        #         query_structure = QUERY_STRUCTURE::PATH;
        #     } else if (degree_array.size() >= 4 && degree_array[0]==1 && degree_array[degree_array.size()-2]==1 && degree_array[degree_array.size()-1]>=3){ // Check if it is a star query...
        #         query_structure = QUERY_STRUCTURE::STAR;
        #     } else if (degree_array.size() >= 6){ // Check if it is a snowflake query...
        #         vector<int>::iterator low_it = lower_bound(degree_array.begin(), degree_array.end(), 3);
        #         if ((degree_array.end()-low_it)>=2 && *(low_it-1)==1){
        #             query_structure = QUERY_STRUCTURE::SNOWFLAKE;
        #         } else {
        #             query_structure = QUERY_STRUCTURE::COMPLEX;
        #         }
        #     } else if (degree_array.size() >= 3){ // Check if it is a complex query...
        #         query_structure = QUERY_STRUCTURE::COMPLEX;
        #     }
        
        #     if (max_join_selectivity >= 0.25){
        #         query_category = QUERY_CATEGORY::HIGH_SELECTIVITY;
        #     } else if (max_join_selectivity < 0.25 && max_join_selectivity >= 0.1){
        #         query_category = QUERY_CATEGORY::MEDIUM_SELECTIVITY;
        #     } else if (max_join_selectivity < 0.1 && max_join_selectivity >= 0.025){
        #         query_category = QUERY_CATEGORY::LOW_SELECTIVITY;
        #     }
        #     */
        """
#     /*
#     pair<QUERY_STRUCTURE::enum_t, QUERY_CATEGORY::enum_t> index_key (query_structure, query_category);
#     if (query_map.find(index_key)==query_map.end()){
#         query_map.insert(pair<pair<QUERY_STRUCTURE::enum_t, QUERY_CATEGORY::enum_t>, map<string, string> >(index_key, map<string, string>()));
#     }
#     if (query_map[index_key].find(query_str)==query_map[index_key].end()){
#         string description = "";
#         description.append("..... ..... ..... ..... ..... ..... ..... ..... ..... .....\n");
#         for (map<string, set<string> >::iterator itr1=variable_map.begin(); itr1!=variable_map.end(); itr1++){
#             description.append(itr1->first);
#             description.append(": ");
#             for (set<string>::iterator itr2=itr1->second.begin(); itr2!=itr1->second.end(); itr2++){
#                 description.append(*itr2);
#                 description.append(" ");
#             }
#             description.append("\n");
#         }
#         description.append("[join-selectivity]:\t");
#         description.append("[");
#         description.append(boost::lexical_cast<string>(min_join_selectivity));
#         description.append(", ");
#         description.append(boost::lexical_cast<string>(max_join_selectivity));
#         description.append("]");
#         description.append("\n");
#         description.append("{ ");
#         for (int i=0; i<degree_array.size(); i++){
#             description.append(boost::lexical_cast<string>(degree_array[i]));
#             description.append(" ");
#         }
#         description.append("}\n");
#         stringstream ss;
#         ss<<query_structure<<"-"<<query_category<<"\n";
#         description.append(ss.str());
#         //description.append(query_structure);
#         //description.append("-");
#         //description.append(query_category);
#         //description.append("\n");
#         description.append("..... ..... ..... ..... ..... ..... ..... ..... ..... .....\n");
#         query_map[index_key].insert(pair<string, string>(query_str, description));
#     */
        pass

#         /// Randomly choose variables, ignore types date, integer, name or string...
#         /// You should also ignore variables corresponding to join vertices...
#         vector<string> eligible_list;
        eligible_list = []
#         for (map<string, set<string> >::iterator itr=variable_map.begin(); itr!=variable_map.end(); itr++){
        for itr_key, itr_value in variable_map.items():
#             string var_type = *(itr->second.begin());
            var_type = itr_value[0]  # #################################################
#             if (var_type.compare("date")==0 || var_type.compare("integer")==0 || var_type.compare("name")==0 || var_type.compare("string")==0){
            if var_type == 'date' or var_type == 'integer' or var_type == 'name' or var_type == 'string':
#                 continue;
                continue
#             }
#             eligible_list.push_back(itr->first);
            eligible_list.append(itr_key)
#         }
#         random_shuffle(eligible_list.begin(), eligible_list.end());
        random.shuffle(eligible_list)

#         /// Select <const-count> variables...
#         if (eligible_list.size()>const_count){
        if len(eligible_list) > const_count:
#             eligible_list.erase(eligible_list.begin()+const_count, eligible_list.end());
            eligible_list = eligible_list[0: const_count-1]
#         }

#         bool varValid = true;
        varValid = True
#         string mapping = "";
        mapping = ""
#         for (int vid=0; vid<eligible_list.size(); vid++){
        for eligible in eligible_list:
#             string var_name = eligible_list[vid];
            var_name = eligible
#             string var_root = var_name.substr(1);
            var_root = var_name[1:]

#             mapping.append("#mapping");
            mapping += '#mapping'
#             mapping.append(" ");
            mapping += ' '
#             mapping.append(var_root); /// Drop preceding '?'...
            mapping += var_root
#             mapping.append(" ");
            mapping += ' '
#             mapping.append(*(variable_map[var_name].begin()));
            mapping += variable_map[var_name][0]
#             mapping.append(" ");
            mapping += ' '
#             mapping.append("uniform");
            mapping += 'uniform'
#             mapping.append("\n");
            mapping += '\n'

#             boost::algorithm::replace_all(query_str, var_name + string(" "), string("%") + var_root + string("%") + string(" "));
            query_str = query_str.replace(f'{var_name} ', f'%{var_root}% ')
#             boost::algorithm::replace_all(query_str, var_name + string("\t"), string("%") + var_root + string("%") + string("\t"));
            query_str = query_str.replace(f'{var_name}\t', f'%{var_root}%\t')

#             if (!constJoinVertexAllowed){
            if not constJoinVertexAllowed:
#                 string placeholder = string("%") + var_root + string("%");
                placeholder = f'%{var_root}%'
#                 if (query_str.find(placeholder)!=string::npos && query_str.find(placeholder)!=query_str.rfind(placeholder)){
                if query_str.find(placeholder) >= 0 and query_str.find(placeholder) != query_str.rfind(placeholder):
#                     varValid = false;
                    varValid = False
#                 }
#             }
#         }
        pass  # end of for

#         bool qValid = false;
        qValid = False
#         string qTemplate = "";
        qTemplate = ''
#         qTemplate.append(mapping);
        qTemplate += mapping
#         qTemplate.append("SELECT ");
        qTemplate += 'SELECT '
#         for (set<string>::iterator itr=variable_set.begin(); itr!=variable_set.end(); itr++){
        for itr in variable_set:
#             string var_name = *itr;
            var_name = itr
#             if (find(eligible_list.begin(), eligible_list.end(), var_name)==eligible_list.end()){
            if var_name not in eligible_list:
#                 qTemplate.append(var_name);
                qTemplate += var_name
#                 qTemplate.append(" ");
                qTemplate += ' '
#                 qValid = true;
                qValid = True
#             }
#         }
#         qTemplate.append("WHERE {");
        qTemplate += 'WHERE {'
#         qTemplate.append("\n");
        qTemplate += '\n'
#         qTemplate.append(query_str);
        qTemplate += query_str
#         qTemplate.append("}");
        qTemplate += '}'
#         qTemplate.append("\n");
        qTemplate += '\n'
#         qTemplate.append("#end");
        qTemplate += '#end'
#         qTemplate.append("\n");
        qTemplate += '\n'

#         if (varValid && qValid){
        if varValid and qValid:
#             cout << qTemplate;
            print(qTemplate)
#             return true;
            return True
#         }
#     /*
#     }
#     */

#     return false;
        return False
# }
        pass  # end of traverse

# void statistics::print_graph() const{
    def print_graph(self):
        """

        Args:

        Returns:

        """
#     int vertex_counter = 0;
        vertex_counter: int = 0
#     cout<<"digraph rdf {"<<"\n";
        print("digraph rdf {")
#     for (map<string, set<statistics_st> >::const_iterator itr1=graph.begin(); itr1!=graph.end(); itr1++){
        for itr1 in self.graph:
#         const set<statistics_st> edge_list = itr1->second;
#         cout<< itr1->first << " " << "[" << "label=\"" << itr1->first << "\"];" << "\n";
#         for (set<statistics_st>::const_iterator itr2=edge_list.begin(); itr2!=edge_list.end(); itr2++){
#             // You need to handle literals separately...
#             if (itr2->_vertex2.compare("date") !=0 && itr2->_vertex2.compare("integer") !=0 && itr2->_vertex2.compare("name") !=0 && itr2->_vertex2.compare("string")){
#                 cout << itr1->first << " -> " << itr2->_vertex2 << " " << "[" << "label=\"" << itr2->_edge << "\"];" << "\n";
#                 //cout << itr2->_vertex2<< " " << "[" << "label=\"" << itr2->_vertex2 << "\"];" << "\n";
#             } else {
#                 string vertex = "v";
#                 vertex.append(boost::lexical_cast<string>(vertex_counter));
#                 cout << itr1->first << " -> " << vertex << " " << "[" << "label=\"" << itr2->_edge << "\"];" << "\n";
#                 vertex_counter++;
#             }
#         }
#     }
            pass  # end of for
#     cout<<"}"<<"\n";
        print("}")
# }
        pass  # end of statistics::print_graph

# /*
# void statistics::compute(){
#     for (map<string, set<statistics_st> >::const_iterator itr1=graph.begin(); itr1!=graph.end(); itr1++){
#         const set<statistics_st> edge_list = itr1->second;
#         for (set<statistics_st>::const_iterator itr2=edge_list.begin(); itr2!=edge_list.end(); itr2++){
#             bool direction = true;
#             if (itr1->first.find(itr2->_vertex1)!=string::npos){
#                 direction = true;
#             } else if (itr1->first.find(itr2->_vertex2)!=string::npos){
#                 direction = false;
#             } else {
#                 cerr<<"Direction is wrong..."<<"\n";
#                 cerr<<"Vertex "<<itr1->first<<"\n";
#                 cerr<<"Edge-1 "<<itr2->_vertex1<<"\n";
#                 cerr<<"Edge-2 "<<itr2->_vertex2<<"\n";
#                 exit(0);
#             }
#             for (int d=0; d<((int)DISTRIBUTION_TYPES::UNDEFINED); d++){
#                 DISTRIBUTION_TYPES::enum_t distribution = (DISTRIBUTION_TYPES::enum_t) d;
#                 string key = get_key(itr1->first, itr2->_edge, direction, distribution);
#                 //cout<<"Sampling "<<key<<"..."<<"\n";
#                 double mean_pr = 0.0;
#                 double mean_card = 0.0;
#                 int sampling_count = 3;
#                 for (int k=0; k<sampling_count; k++){
#                     pair<double, double> stats = sample(itr1->first, distribution, itr2->_edge, direction);
#                     mean_pr += stats.first;
#                     mean_card += stats.second;
#                 }
#                 mean_pr = mean_pr / ((double) sampling_count);
#                 mean_card = mean_card / ((double) sampling_count);
#                 _statistics_table.insert(pair<string, pair<double, double> >(key, pair<double, double>(mean_pr, mean_card)));
#             }
#         }
#     }
# }
# */

# /*
# void statistics::print_stats() const{
#     for (map<string, pair<double, double> >::const_iterator itr1=_statistics_table.begin(); itr1!=_statistics_table.end(); itr1++){
#         string key = itr1->first;
#         pair<double, double> value = itr1->second;
#         cout<<key<<"\t"<<value.first<<"\t"<<value.second<<"\n";
#     }
# }
# */

# void statistics::index_triples(const vector<triple_st> & triple_array){
    def index_triples(self, triple_array):
        """

        Args:

        Returns:

        """
#     for (int i=0; i<triple_array.size(); i++){
        for triple in triple_array:
#         _spo_index.push_back(triple_array[i]);
            self._spo_index.append(triple)
#         _ops_index.push_back(triple_array[i]);
            self._ops_index.append(triple)
#     }
            pass  # end of for
#     sort(_spo_index.begin(), _spo_index.end(), s_compare());
        self._spo_index.sort(key=lambda x: x._subject)
#     sort(_ops_index.begin(), _ops_index.end(), o_compare());
        self._ops_index.sort(key=lambda x: x._object)
# }
        pass  # end of statistics::index_triples

# string statistics::get_key (string entity, string predicate, bool direction, DISTRIBUTION_TYPES::enum_t distribution) const{
    def get_key(self, entity, predicate, direction, distribution):
        """

        Args:

        Returns:

        """
#     string result = "";
        result: str = ""
#     result.append(entity);
        result += entity
#     result.append("#");
        result += "#"
#     result.append(predicate);
        result += predicate
#     result.append("#");
        result += "#"
#     if (direction){
        if direction:
#         result.append("t");
            result += "t"
#     } else {
        else:
#         result.append("f");
            result += "f"
#     }
        pass  # end of if
#     result.append("#");
        result += "#"
#     switch (distribution){
#         case DISTRIBUTION_TYPES::UNIFORM:{
        if distribution == DistributionTypes.UNIFORM:
#             result.append("uniform");
            result += "uniform"
#             break;
#         }
#         case DISTRIBUTION_TYPES::NORMAL:{
        elif distribution == DistributionTypes.NORMAL:
#             result.append("normal");
            result += "normal"
#             break;
#         }
#         case DISTRIBUTION_TYPES::ZIPFIAN:{
        elif distribution == DistributionTypes.ZIPFIAN:
#             result.append("zipfian");
            result += "zipfian"
#             break;
#         }
#         case DISTRIBUTION_TYPES::UNDEFINED:{
        elif distribution == DistributionTypes.UNDEFINED:
#             result.append("undefined");
            result += "undefined"
#             break;
#         }
        pass  # end of if
#     }
#     return result;
        return result
# }
        pass  # end of statistics::get_key

        # /*
        # pair<double, double> statistics::sample (string entity, DISTRIBUTION_TYPES::enum_t distribution, string predicate, bool direction) const{
        #     mapping_m_t * generator = NULL;
        #     int pos = string::npos;
        #     if ((pos=entity.find("@"))!=string::npos){
        #         string type = entity.substr(0, pos);
        #         string restriction = entity.substr(pos+1);
        #         generator = new mapping_m_t (string("?x"), type, restriction, distribution);
        #     } else {
        #         generator = new mapping_m_t (string("?x"), entity, distribution);
        #     }
        #     double pr_exists = 0.0;
        #     double mean_card = 0.0;
        #     unsigned int instance_count = 0;
        #     generator->generate(*_model, instance_count);
        #     int sampling_factor = instance_count * 5;
        #     for (int i=0; i<sampling_factor; i++){
        #         int cardinality = 0;
        #         string min_subject = "", min_predicate = "", min_object = "";
        #         min_predicate.append("<");
        #         min_predicate.append(_model->_namespace_map.replace(predicate));
        #         min_predicate.append(">");
        #         string rdf_term = generator->generate(*_model);
        #         if (direction){
        #             min_subject.append(rdf_term);
        #             triple_st min_range (min_subject, min_predicate, min_object);
        #             //cout<<"Searching for "<<min_range<<"\n";
        #             vector<triple_st>::const_iterator itr = lower_bound(_spo_index.begin(), _spo_index.end(), min_range, s_compare());
        #             for (; itr!=_spo_index.end() && itr->_subject.compare(min_subject)==0; itr++){
        #                 //cout<<"Debug "<<*itr<<"\n";
        #                 if (itr->_predicate.compare(min_predicate)==0){
        #                     cardinality++;
        #                 }
        #             }
        #         } else {
        #             min_object.append(rdf_term);
        #             triple_st min_range (min_subject, min_predicate, min_object);
        #             //cout<<"Searching for "<<min_range<<"\n";
        #             vector<triple_st>::const_iterator itr = lower_bound(_ops_index.begin(), _ops_index.end(), min_range, o_compare());
        #             for (; itr!=_ops_index.end() && itr->_object.compare(min_object)==0; itr++){
        #                 //cout<<"Debug "<<*itr<<"\n";
        #                 if (itr->_predicate.compare(min_predicate)==0){
        #                     cardinality++;
        #                 }
        #             }
        #         }
        #         if (cardinality>0){
        #             pr_exists += 1.0;
        #         }
        #         mean_card += (double) cardinality;
        #     }
        #     mean_card = mean_card / (pr_exists+0.0000000001);
        #     pr_exists = pr_exists / ((double) sampling_factor);
        #     pair<double, double> result = pair<double, double> (pr_exists, mean_card);
        #     return result;
        #     //cout << "\t" << "[" << pr_exists << ", " << mean_card << "]" << "\n";
        # }
        # */

    pass  # end of class Statistics

# -------------------------------------------------------------------------------------------------------

# static unsigned int MAX_LOOP_COUNTER = 50;
MAX_LOOP_COUNTER = 50
# static int MAX_LITERAL_WORDS = 25;
MAX_LITERAL_WORDS = 25

# static map<int,vector<double>*> zipfian_cache;
zipfian_cache = {}

# static boost::mt19937 BOOST_RND_GEN = boost::mt19937(static_cast<unsigned> (time(0)));
# static boost::normal_distribution<double> BOOST_NORMAL_DIST = boost::normal_distribution<double>(0.5, (0.5/3.0));
# static boost::variate_generator<boost::mt19937, boost::normal_distribution<double> > BOOST_NORMAL_DIST_GEN (BOOST_RND_GEN, BOOST_NORMAL_DIST);

# using namespace std;

# ostream& operator<<(ostream& os, const DISTRIBUTION_TYPES::enum_t & distribution){
#     switch (distribution){
#         case DISTRIBUTION_TYPES::UNIFORM:{
#             os << "uniform";
#             break;
#         }
#         case DISTRIBUTION_TYPES::NORMAL:{
#             os << "normal";
#             break;
#         }
#         case DISTRIBUTION_TYPES::ZIPFIAN:{
#             os << "zipfian";
#             break;
#         }
#         case DISTRIBUTION_TYPES::UNDEFINED:{
#             os << "undefined";
#             break;
#         }
#     }
#     return os;
# }

# namespace LITERAL_TYPES {
#     enum enum_t {
#         INTEGER,
#         STRING,
#         NAME,
#         DATE,
#         UNDEFINED
#     };
# };


class LiteralTypes(Enum):
    INTEGER = auto
    STRING = auto
    NAME = auto
    DATE = auto
    UNDEFINED = auto

# namespace DISTRIBUTION_TYPES {
#     enum enum_t {
#         UNIFORM,
#         NORMAL,
#         ZIPFIAN,
#         DYNAMIC,
#         UNDEFINED
#     };
# };


class DistributionTypes(Enum):
    UNIFORM = 1
    NORMAL = 2
    ZIPFIAN = 3
    DYNAMIC = 4
    UNDEFINED = 5

# ostream& operator<<(ostream& os, const DISTRIBUTION_TYPES::enum_t & distribution);

# namespace OPERATION_TYPES {
#     enum enum_t {
#         ADDITION,
#         MULTIPLICATION,
#         MOD,
#         UNDEFINED
#     };
# };


class OperationTypes(Enum):
    ADDITION = auto
    MULTIPLICATION = auto
    MOD = auto
    UNDEFINED = auto

# ostream& operator<<(ostream& os, const triple_st & triple);

# struct s_compare {
#     bool operator() (const triple_st & lhs, const triple_st & rhs) const;
# };

# struct o_compare {
#     bool operator() (const triple_st & lhs, const triple_st & rhs) const;
# };


def parse_line(line: str) -> list[str]:
    """

    Args:

    Returns:

    """
    xxx = line.replace('\t\n', '')
    xxx = xxx.replace('\n', '')
    xxx = xxx.replace('\t', ' ')
    xxx = xxx.replace('  ', ' ')
    xxx = xxx.replace('  ', ' ')
    xxx = xxx.replace('  ', ' ')
    xxx.strip()
    return xxx.split(' ')


# struct triple_st {
# };
class TripleST:
    """
    Stores a triple

    Attributes:
        _subject (str): subject of a triple
        _predicate (str): predicate of a triple
        _object (str): object of a triple

    """
# triple_st::triple_st (const string & line){
    def __init__(self, line=None, subject=None, predicate=None, object0=None):
        """

        Args:

        Returns:

        """
#     string _subject;
        self._subject: str = ''
#     string _predicate;
        self._predicate: str = ''
#     string _object;
        self._object: str = ''
#
#     triple_st (const string & line);
#     triple_st (const string & subject, const string & predicate, const string & object);
#     bool operator== (const triple_st & rhs) const;
#
#     static vector<triple_st> parse_file (const char * filename);
#     vector<string> result;
        self.result = []
        if line is not None:
#     string trimmed = line.substr(0, line.find_last_of("."));
            trimmed = line[:line.rfind('.')]
#     boost::trim(trimmed);
            trimmed = trimmed.strip()
#     boost::algorithm::split(result, trimmed, boost::is_any_of("\t"));
            result = trimmed.split('\t')
#     _subject = result[0];
            self._subject = result[0]
#     _predicate = result[1];
            self._predicate = result[1]
#     _object = result[2];
            self._object = result[2]
# }

# triple_st::triple_st (const string & subject, const string & predicate, const string & object){
        if subject is not None and predicate is not None and object0 is not None:
#     _subject = subject;
            self._subject = subject
#     _predicate = predicate;
            self._predicate = predicate
#     _object = object;
            self._object = object0
# }

# bool triple_st::operator== (const triple_st & rhs) const{
#     return _subject.compare(rhs._subject)==0 && _predicate.compare(rhs._predicate)==0 && _object.compare(rhs._object)==0;
# }

    @staticmethod
# vector<triple_st> triple_st::parse_file (const char * filename){
    def parse_file(filename):
        """
        Read triples from an RDF file

        Args:
            filename (str): path to an RDF file
        Returns:

        """
#     vector<triple_st> result;
#     ifstream ifs (filename);
        with open(filename, 'r') as ifs:  # open an RDF file
            lines = ifs.readlines()  # read all the lines in the file
            result: list[str] = [None] * len(lines)  # prepare a list
#     string line;
#     while ( getline(ifs, line) ){
            for index, line in enumerate(lines):
#         triple_st cur_triple (line);
                cur_triple: TripleST = TripleST(line=line)
#         result.push_back(cur_triple);
                result[index] = cur_triple
#         cout << cur_triple << "\n";
#                 print(f'{cur_triple._subject}\t{cur_triple._predicate}\t{cur_triple._object}')
                if index % 100000 == 0:
                    print('.', end='')  # progress indicator, because this step takes long time

#     }
            print()  # CR
#     ifs.close();
#     return result;
        return result
# }
        pass  # end of parse_file
    pass  # end of TripleST


# ostream& operator<<(ostream& os, const triple_st & triple){
#     os<<triple._subject<<"\t"<<triple._predicate<<"\t"<<triple._object;
#     return os;
# }

# bool s_compare::operator() (const triple_st & lhs, const triple_st & rhs) const{
#     if (lhs._subject.compare(rhs._subject)!=0){
#         return lhs._subject.compare(rhs._subject) < 0;
#     }
#     if (lhs._predicate.compare(rhs._predicate)!=0){
#         return lhs._predicate.compare(rhs._predicate) < 0;
#     }
#     if (lhs._object.compare(rhs._object)!=0){
#         return lhs._object.compare(rhs._object) < 0;
#     }
#     return false;
# }

# bool o_compare::operator() (const triple_st & lhs, const triple_st & rhs) const{
#     if (lhs._object.compare(rhs._object)!=0){
#         return lhs._object.compare(rhs._object) < 0;
#     }
#     if (lhs._predicate.compare(rhs._predicate)!=0){
#         return lhs._predicate.compare(rhs._predicate) < 0;
#     }
#     if (lhs._subject.compare(rhs._subject)!=0){
#         return lhs._subject.compare(rhs._subject) < 0;
#     }
#     return false;
# }

# struct namespace_m_t {
class NamespaceMT:
    """

    Attributes:
        _alias (str):
        _prefix (str):
    """

#     namespace_m_t (string token);

#     static namespace_m_t * parse (const string & line);
# };
# namespace_m_t::namespace_m_t (string token){
    def __init__(self, token: str):
        """

        Args:

        Returns:

        """
#     string                          _alias;
        self._alias = ''
#     string                          _prefix;
        self._prefix = ''
#     unsigned int pos = token.find_first_of('=');
        pos = token.find('=')
        token_split = token.split('=')
#     if (pos!=string::npos){
        if pos >= 0:
#         _alias = token.substr(0, pos);
            self._alias = token_split[0]
#         _prefix = token.substr(pos+1);
            self._prefix = token_split[1]
#         //cout<<"namespace ["<<_alias<<"]-->["<<_prefix<<"] generated..."<<"\n";
#     } else {
        else:
#         cerr<<"[namespace_m_t::namespace_m_t] Expected [alias]:[prefix]..."<<"\n";
            print('[namespace_m_t::namespace_m_t] Expected [alias]:[prefix]...')
#         exit(0);
            sys.exit(0)
#     }
# }
        pass  # end of namespace_m_t::namespace_m_t

    @staticmethod
# namespace_m_t * namespace_m_t::parse (const string & line){
    def parse(line):
        """

        Args:

        Returns:

        """
#     string namespace_declaration;
        namespace_declaration: str = ''
#
#     stringstream parser(line);
        parser = parse_line(line)
#     int index = 0;
        index: int = 0
#     string token;
#
#     while (parser>>token){
        for token in parser:
#         switch (index){
#             case 0:{
            if index == 0:
#                 if (token.compare("#namespace")!=0){
                if token.find('#namespace') < 0:
#                     cerr<<"[predicate_m_t::parse()]\tExpecting #namespace..."<<"\n";
                    print('[predicate_m_t::parse()]\tExpecting #namespace...')
#                     exit(0);
                    sys.exit(0)
#                 }
#                 break;
#             }
#             case 1:{
            if index == 1:
#                 namespace_declaration = token;
                namespace_declaration = token
#                 break;
#             }
#         }
#         index++;
            index += 1
#     }
#
#     if (index==2){
        if index == 2:
#         return new namespace_m_t(namespace_declaration);
            return NamespaceMT(namespace_declaration)
#     } else {
        else:
#         cerr<<"[namespace_m_t::parse()]\tExpecting 1 argument..."<<"\n";
            print('[namespace_m_t::parse()]\tExpecting 1 argument...')
#         exit(0);
            sys.exit(0)
#     }
# }
        pass  # end of

# void namespace_map::to_str (vector<string> & lines) const{
    def to_str(self, lines):
        """

        Args:

        Returns:

        """
#     for (map<string, string>::const_iterator itr=_index.begin(); itr!=_index.end(); itr++){
        for itr in self._index:
#         string line = "";
#         line.append(itr->first);
#         line.append(" ");
#         line.append(itr->second);
            line = f'{itr.first} {itr.second}'
#         lines.push_back(line);
            lines.append(line)
#     }
# }
        pass  # end of to_str
    pass  # end of NamespaceMT


class TypeMap:
    """

    Attributes:

    """
# class type_map {
#     public:
#         type_map();
#         ~type_map();

#         void clear();
#         void insert (const string & instance, const string & type);
#         bool instanceof (const string & instance, const string & type) const;
#         vector<string> * get_instances (const string & entity, const string & type) const;

#         void print () const;
#         void to_str (vector<string> & lines) const;
#     private:
#         /// map<string, set<string> > _index;
# };
# type_map::type_map(){
    def __init__(self):
        """

        Args:

        Returns:

        """
#         unordered_map<string, unordered_set<string> > _index; // @suppress("Invalid template argument")
        self._index = {}
#
# }

# type_map::~type_map(){
#
# }

# void type_map::clear(){
    def clear(self):
        """

        Args:

        Returns:

        """
    #     _index.clear();
        self._index = {}
# }

# void type_map::insert (const string & instance, const string & type){
    def insert(self, instance, type0):
        """

        Args:

        Returns:

        """
#     if (_index.find(type)==_index.end()){
        try:
            xxx = self._index[type0]
#         _index.insert(pair<string, unordered_set<string> >(type, unordered_set<string>()));
        except KeyError:
            self._index[type0] = set()
#     }
#     _index[type].insert(instance);
        self._index[type0].add(instance)
# }
        pass  # end of clear

# bool type_map::instanceof (const string & instance, const string & type) const{
    def instance_of(self, instance, type0):
        """

        Args:

        Returns:

        """
#     unordered_map<string, unordered_set<string> >::const_iterator f_it1 = _index.find(type);
        try:
            f_it1 = self._index[type0]
#     if (f_it1!=_index.end()){
#         unordered_set<string>::const_iterator f_it2 = f_it1->second.find(instance);
#         if (f_it2!=f_it1->second.end()){
            if instance in f_it1:
#             return true;
                return True
#         }
#     }
        except KeyError:
#     return false;
            pass
        return False
# }
        pass  # end of instance_of

# vector<string> * type_map::get_instances (const string & entity, const string & type) const{
    def get_instances(self, entity, type0):
        """

        Args:

        Returns:

        """
#     //cout<<"Looking up type restriction:"<<type<<"\n";
#     unordered_map<string, unordered_set<string> >::const_iterator f_it1=_index.find(type);
        try:
            f_it1 = self._index[type0]
#     if (f_it1!=_index.end()){
#         vector<string> * result = new vector<string>();
#         for (unordered_set<string>::const_iterator f_it2=f_it1->second.begin(); f_it2!=f_it1->second.end(); f_it2++){
#             string instance = *f_it2;
#             if (boost::starts_with(instance, entity)){
#                 result->push_back(*f_it2);
#             }
#         }
            result = list[f_it1]
#         return result;
            return result
#     }
        except KeyError:
            pass
#     //cout<<"Type restriction not found..."<<"\n";
#     //cout<<"Listing contents under type_map..."<<"\n";
#     //for (map<string, set<string> >::const_iterator itr1=_index.begin(); itr1!=_index.end(); itr1++){
#         //cout<<"\t"<<itr1->first<<"\n";
#     //}
#     //cout<<"Contents listed..."<<"\n";
#     return NULL;
        return None
# }
        pass  # end of get_instances

# void type_map::print() const{
    def print(self):
        """

        Args:

        Returns:

        """
#     for (unordered_map<string, unordered_set<string> >::const_iterator itr1=_index.begin(); itr1!=_index.end(); itr1++){
        for itr1_key, itr1_value in self._index.items():
#         cout<<"Type:::\t"<<itr1->first<<"\n";
            print(f'Type:::\t{itr1_key}')
#         for (unordered_set<string>::const_iterator itr2=itr1->second.begin(); itr2!=itr1->second.end(); itr2++){
            for itr2 in itr1_value:
#             cout<<"\t\t\t"<<*itr2<<"\n";
                print(f'\t\t\t{itr2}')
#         }
#     }
# }
        pass  # end of print

# void type_map::to_str (vector<string> & lines) const{
    def to_str(self):
        """

        Args:

        Returns:

        """
        lines = []
#     for (unordered_map<string, unordered_set<string> >::const_iterator itr1=_index.begin(); itr1!=_index.end(); itr1++){
        for itr1_key, itr1_value in self._index.items():
#         string line = "";
            line: str = ''
#         line.append(itr1->first);
            line += itr1_key
#         line.append(" ");
            line += ' '
#         for (unordered_set<string>::const_iterator itr2=itr1->second.begin(); itr2!=itr1->second.end(); itr2++){
            for itr2 in itr1_value:
#             line.append(*itr2);
                line += itr2
#             line.append(" ");
                line += ' '
#         }
#         lines.push_back(line);
            lines.append(line)
#     }
        return lines
# }
        pass  # end of type_map::to_str
    pass  # end of TypeMap


# struct predicate_m_t {
class PredicateMT:
    """

    Attributes:

    """
#     string                          _label;
#     LITERAL_TYPES::enum_t           _literal_type;
#     string                          _range_min;
#     string                          _range_max;
#     DISTRIBUTION_TYPES::enum_t      _distribution_type;

#     void init (string label, LITERAL_TYPES::enum_t literal_type);

#     predicate_m_t (string label, LITERAL_TYPES::enum_t literal_type);
#     predicate_m_t (string label, LITERAL_TYPES::enum_t literal_type, string range_min, string range_max);
#     predicate_m_t (string label, LITERAL_TYPES::enum_t literal_type, string range_min, string range_max, DISTRIBUTION_TYPES::enum_t distribution_type);
#     predicate_m_t (const predicate_m_t & rhs);

#     string generate (const namespace_map & n_map);

#     static predicate_m_t * parse (const string & line);
# };

# void predicate_m_t::init(string label, LITERAL_TYPES::enum_t literal_type){
    def init(self, label, literal_type):
        """

        Args:

        Returns:

        """
#     _label = label;
        self._label = label
#     _literal_type = literal_type;
        self._literal_type = literal_type
#     switch (_literal_type){
#         case LITERAL_TYPES::INTEGER:{
        if self._literal_type == LiteralTypes.INTEGER:
#             _range_min = boost::lexical_cast<string>(numeric_limits<unsigned short>::min());
            self._range_min = str(-1000000)
#             _range_max = boost::lexical_cast<string>(numeric_limits<unsigned short>::max());
            self._range_max = str(1000000)
#             break;
#         }
#         case LITERAL_TYPES::STRING:
#         case LITERAL_TYPES::NAME:{
        elif self._literal_type == LiteralTypes.STRING or self._literal_type == LiteralTypes.NAME:
#             _range_min = string("A");
            self._range_min = 'A'
#             _range_max = string("z");
            self._range_max = 'z'
#             break;
#         }
#         case LITERAL_TYPES::DATE:{
        elif self._literal_type == LiteralTypes.DATE:
#             boost::posix_time::ptime cur_time (boost::posix_time::second_clock::local_time());
            cur_time = datetime.now().time()
#             boost::gregorian::date cur_date = cur_time.date();
            cur_date = datetime.now().date()
#             _range_min = string("1970-01-01");
            self._range_min = str('1970-01-01')
#             _range_max = boost::gregorian::to_iso_extended_string(cur_date);
            self._range_max = str(cur_date)
#             break;
#         }
#     }
#     _distribution_type = DISTRIBUTION_TYPES::UNIFORM;
        self._distribution_type = DistributionTypes.UNIFORM
# }

# predicate_m_t::predicate_m_t (string label, LITERAL_TYPES::enum_t literal_type){
    def __init__(self, label='', literal_type=None, range_min='', range_max='', distribution_type=None, rhs=None):
        """

        Args:

        Returns:

        """
#     string                          _label;
        self._label: str = ''
#     LITERAL_TYPES::enum_t           _literal_type;
        self._literal_type = None
#     string                          _range_min;
        self._range_min: str = ''
#     string                          _range_max;
        self._range_max: str = ''
#     DISTRIBUTION_TYPES::enum_t      _distribution_type;
        self._distribution_type = None
#     init(label, literal_type);
        if label != '' and literal_type is not None:
            self.init(label, literal_type)
# }
        pass  # end of predicate_m_t::predicate_m_t (string label, LITERAL_TYPES::enum_t literal_type)

# predicate_m_t::predicate_m_t (string label, LITERAL_TYPES::enum_t literal_type, string range_min, string range_max){
#     init(label, literal_type);
#     _range_min = range_min;
        if range_min != '':
            self._range_min = range_min
#     _range_max = range_max;
        if range_max != '':
            self._range_max = range_max
# }
        pass  # end of predicate_m_t::predicate_m_t (string label, LITERAL_TYPES::enum_t literal_type, string range_min, string range_max)

# predicate_m_t::predicate_m_t (string label, LITERAL_TYPES::enum_t literal_type, string range_min, string range_max, DISTRIBUTION_TYPES::enum_t distribution_type){
#     init(label, literal_type);
#     _range_min = range_min;
#     _range_max = range_max;
#     _distribution_type = distribution_type;
        if distribution_type is not None:
            self._distribution_type = distribution_type
# }
        pass  # end of predicate_m_t::predicate_m_t (string label, LITERAL_TYPES::enum_t literal_type, string range_min, string range_max, DISTRIBUTION_TYPES::enum_t distribution_type)

# predicate_m_t::predicate_m_t (const predicate_m_t & rhs){
        if rhs is not None:
#     _label = rhs._label;
            self._label = rhs._label
#     _literal_type = rhs._literal_type;
            self._literal_type = rhs._literal_type
#     _range_min = rhs._range_min;
            self._range_min = rhs._range_min
#     _range_max = rhs._range_max;
            self._range_max = rhs._range_max
#     _distribution_type = rhs._distribution_type;
            self._distribution_type = rhs._distribution_type
# }
        pass  # end of predicate_m_t::predicate_m_t (const predicate_m_t & rhs)

    @staticmethod
# predicate_m_t * predicate_m_t::parse (const string & line){
    def parse(line):
        """

        Args:

        Returns:

        """
#     string label;
        label = ''
#     LITERAL_TYPES::enum_t literal_type = LITERAL_TYPES::UNDEFINED;
        literal_type = LiteralTypes.UNDEFINED
#     string range_min;
        range_min = ''
#     string range_max;
        range_max = ''
#     DISTRIBUTION_TYPES::enum_t distribution_type = DISTRIBUTION_TYPES::UNDEFINED;
        distribution_type = DistributionTypes.UNDEFINED
#
#     stringstream parser(line);
        parser = parse_line(line)
#     int index = 0;
        index: int = 0
#     string token;
#
#     while (parser>>token){
        for token in parser:
#         switch (index){
#             case 0:{
            if index == 0:
#                 if (token.compare("#predicate")!=0){
                if token != '#predicate':
#                     cerr<<"[predicate_m_t::parse()]\tExpecting #predicate..."<<"\n";
                    print('[predicate_m_t::parse()]\tExpecting #predicate...')
#                     exit(0);
                    sys.exit(0)
#                 }
#                 break;
#             }
#             case 1:{
            elif index == 1:
#                 label = token;
                label = token
#                 break;
#             }
#             case 2:{
            elif index == 2:
#                 if (token.compare("integer")==0 || token.compare("INTEGER")==0){
                if token == 'integer' or token == 'INTEGER':
#                     literal_type = LITERAL_TYPES::INTEGER;
                    literal_type = LiteralTypes.INTEGER
#                 } else if (token.compare("string")==0 || token.compare("STRING")==0){
                if token == 'string' or token == 'STRING':
#                     literal_type = LITERAL_TYPES::STRING;
                    literal_type = LiteralTypes.STRING
#                 } else if (token.compare("name")==0 || token.compare("NAME")==0){
                if token == 'name' or token == 'NAME':
#                     literal_type = LITERAL_TYPES::NAME;
                    literal_type = LiteralTypes.NAME
#                 } else if (token.compare("date")==0 || token.compare("DATE")==0){
                if token == 'date' or token == 'DATE':
#                     literal_type = LITERAL_TYPES::DATE;
                    literal_type = LiteralTypes.DATE
#                 }
#                 break;
#             }
#             case 3:{
            elif index == 3:
#                 range_min = token;
                range_min = token
#                 break;
#             }
#             case 4:{
            elif index == 4:
#                 range_max = token;
                range_max = token
#                 break;
#             }
#             case 5:{
            elif index == 5:
#                 if (token.compare("uniform")==0 || token.compare("UNIFORM")==0){
                if token == 'uniform' or token == 'UNIFORM':
#                     distribution_type = DISTRIBUTION_TYPES::UNIFORM;
                    distribution_type = DistributionTypes.UNIFORM
#                 } else if (token.compare("normal")==0 || token.compare("NORMAL")==0){
                if token == 'normal' or token == 'NORMAL':
#                     distribution_type = DISTRIBUTION_TYPES::NORMAL;
                    distribution_type = DistributionTypes.NORMAL
#                 } else if (token.compare("zipfian")==0 || token.compare("ZIPFIAN")==0){
                if token == 'zipfian' or token == 'ZIPFIAN':
#                     distribution_type = DISTRIBUTION_TYPES::ZIPFIAN;
                    distribution_type = DistributionTypes.ZIPFIAN
#                 }
#                 break;
                break
#             }
#         }
#         index++;
            index += 1
#     }

#     if (index==3){
        if index == 3:
#         return new predicate_m_t(label, literal_type);
            return PredicateMT(label=label, literal_type=literal_type)
#     } else if (index==5){
        if index == 5:
#         return new predicate_m_t(label, literal_type, range_min, range_max);
            return PredicateMT(label=label, literal_type=literal_type, range_min=range_min, range_max=range_max)
#     } else if (index==6){
        if index == 5:
#         return new predicate_m_t(label, literal_type, range_min, range_max, distribution_type);
            return PredicateMT(label=label, literal_type=literal_type, range_min=range_min, range_max=range_max, distribution_type=distribution_type)
#     } else {
        else:
#         cerr<<"[predicate_m_t::parse()]\tExpecting 2, 3 or 5 arguments..."<<"\n";
            print('[predicate_m_t::parse()]\tExpecting 2, 3 or 5 arguments...')
#         exit(0);
            sys.exit(0)
#     }
        pass  # end of if
# }
        pass  # end of parse

# string predicate_m_t::generate (const namespace_map & n_map){
    def generate(self, n_map):
        """

        Args:

        Returns:

        """
#     string result = "";
        result: str = ''
#     string literal = model::generate_literal(_literal_type, _distribution_type, _range_min, _range_max);
        literal = Model.generate_literal(self._literal_type, self._distribution_type, self._range_min, self._range_max)
#     result.append("<");
        result += '<'
#     result.append(n_map.replace(_label));
        result += n_map.replace(self._label)
#     result.append(">");
        result += '>'
#     result.append("\t");
        result += '\t'
#     result.append("\"");
        result += '"'
#     result.append(literal);
        result += literal
#     result.append("\"");
        result += '"'
#     return result;
        return result
# }
        pass  # end of generate
    pass  # end of PredicateMT


# struct predicate_group_m_t {
class PredicateGroupMT:
    """

    Attributes:

    """

    #     predicate_group_m_t ();
    #     predicate_group_m_t (float gen_probability);
    #     predicate_group_m_t (float gen_probability, const string & type_restriction);
    #     predicate_group_m_t (const predicate_group_m_t & rhs);
    #     ~predicate_group_m_t();

    #     static predicate_group_m_t * parse (const string & line);
    # };
    # predicate_group_m_t::predicate_group_m_t (){
#     _post_process = false;
#     _gen_probability = 1.0;
#     _type_restriction = NULL;
# }

# predicate_group_m_t::predicate_group_m_t (float gen_probability){
    def __init__(self, gen_probability=None, type_restriction=None, rhs=None):
        """

        Args:

        Returns:

        """
        #     bool                            _post_process;
        self._post_process: bool = False
    #     float                           _gen_probability;
        self._gen_probability: float = 0.0
    #     string *                        _type_restriction;
        self._type_restriction: str = ''
    #     vector<predicate_m_t*>          _predicate_array;
        self._predicate_array = []
    #     _post_process = false;
#     _gen_probability = gen_probability;
        if gen_probability is not None:
            self._gen_probability = gen_probability
#     _type_restriction = NULL;
# }
        pass  # end of predicate_group_m_t::predicate_group_m_t (float gen_probability)

# predicate_group_m_t::predicate_group_m_t (float gen_probability, const string & type_restriction){
        if gen_probability is not None and type_restriction != '':
#     _post_process = true;
            self._post_process = True
#     _gen_probability = gen_probability;
            self._gen_probability = gen_probability
#     _type_restriction = new string(type_restriction);
            self._type_restriction = type_restriction
# }
        pass  # end of predicate_group_m_t::predicate_group_m_t (float gen_probability, const string & type_restriction)

# predicate_group_m_t::predicate_group_m_t (const predicate_group_m_t & rhs){
        if rhs is not None:
#     _post_process = rhs._post_process;
            self._post_process = rhs._post_process
#     _gen_probability = rhs._gen_probability;
            self._gen_probability = rhs._gen_probability
#     _type_restriction = new string(*(rhs._type_restriction));
            self._type_restriction = rhs._type_restriction
#     for (unsigned int i=0; i<rhs._predicate_array.size(); i++){
            for ppp in rhs._predicate_array:
#         _predicate_array.push_back(new predicate_m_t(*(rhs._predicate_array[i])));
                self._predicate_array.append(PredicateMT(rhs=ppp))
#     }
                pass  # end of for
# }
        pass  # end of predicate_group_m_t::predicate_group_m_t (const predicate_group_m_t & rhs)

# predicate_group_m_t::~predicate_group_m_t(){
#     for (unsigned int i=0; i<_predicate_array.size(); i++){
#         delete _predicate_array[i];
#     }
#     delete _type_restriction;
# }

    @staticmethod
# predicate_group_m_t * predicate_group_m_t::parse (const string & line){
    def parse(line):
        """

        Args:

        Returns:

        """
#     float gen_probability = 1.0;
        gen_probability: float = 1.0
#     string type_restriction;
        type_restriction = ''

#     stringstream parser(line);
        parser = parse_line(line)
#     int index = 0;
        index: int = 0
#     string token;
#     while (parser>>token){
        for token in parser:
#         switch (index){
#             case 0:{
            if index == 0:
#                 if (token.compare("<pgroup>")!=0){
                if token != '<pgroup>':
#                     cerr<<"[predicate_group_m_t::parse()]\tExpecting <pgroup>..."<<"\n";
                    print('[predicate_group_m_t::parse()]\tExpecting <pgroup>...')
#                     exit(0);
                    sys.exit(0)
#                 }
#                 break;
#             }
#             case 1:{
            if index == 1:
#                 gen_probability = boost::lexical_cast<float>(token);
                gen_probability = float(token)
#                 break;
#             }
#             case 2:{
            if index == 2:
#                 if (!boost::starts_with(token, "@")){
                if not token.startswith('@'):
#                     cerr<<"[predicate_group_m_t::parse()]\tExpecting '@' qualifier before type restriction..."<<"\n";
                    print('[predicate_group_m_t::parse()]\tExpecting '@' qualifier before type restriction...')
#                     exit(0);
                    sys.exit(0)
#                 }
#                 type_restriction = token.substr(1);
                type_restriction = token[1:]
#                 break;
#             }
#         }
#         index++;
            index += 1
#     }

#     if (index==1){
        if index == 1:
#         return new predicate_group_m_t();
            return PredicateGroupMT()
#     } else if (index==2){
        elif index == 2:
#         return new predicate_group_m_t(gen_probability);
            return PredicateGroupMT(gen_probability=gen_probability)
#     } else if (index==3){
        elif index == 3:
#         return new predicate_group_m_t(gen_probability, type_restriction);
            return PredicateGroupMT(gen_probability=gen_probability, type_restriction=type_restriction)
#     } else {
        else:
#         cerr<<"[predicate_group_m_t::parse()]\tExpecting 0, 1 or 2 arguments..."<<"\n";
            print('[predicate_group_m_t::parse()]\tExpecting 0, 1 or 2 arguments...')
#         exit(0);
            sys.exit(0)
#     }
# }
    pass  # end of PredicateGroupMT


# struct resource_m_t {
class ResourceMT:
    """

    Attributes:

    """
#     resource_m_t (bool scalable, string type_prefix, unsigned int scaling_coefficient);
#     resource_m_t (const resource_m_t & rhs);
#     ~resource_m_t ();

#     void generate (const namespace_map & n_map, map<string, unsigned int> & id_cursor_map); // @suppress("Invalid template argument")
#     void process_type_restrictions (const namespace_map & n_map, const type_map & t_map, const map<string, unsigned int> & id_cursor_map); // @suppress("Invalid template argument")

#     static resource_m_t * parse (const string & line);
# };

# resource_m_t::resource_m_t (bool scalable, string type_prefix, unsigned int scaling_coefficient){
    def __init__(self, scalable=None, type_prefix=None, scaling_coefficient=None, rhs=None):
        """

        Args:

        Returns:

        """
#     bool                            _scalable;
        self._scalable = False
#     string                          _type_prefix;
        self._type_prefix = ''
#     unsigned int                    _scaling_coefficient;
        self._scaling_coefficient = 0
#     vector<predicate_group_m_t*>    _predicate_group_array;
        self._predicate_group_array = []
        if scalable is not None and type_prefix is not None and scaling_coefficient is not None:
#     _scalable = scalable;
            self._scalable = scalable
#     _type_prefix = type_prefix;
            self._type_prefix = type_prefix
#     _scaling_coefficient = scaling_coefficient;
            self._scaling_coefficient = scaling_coefficient
# }
        pass  # end of resource_m_t::resource_m_t (bool scalable, string type_prefix, unsigned int scaling_coefficient)

# resource_m_t::resource_m_t (const resource_m_t & rhs){
        if rhs is not None:
#     _scalable = rhs._scalable;
            self._scalable = rhs._scalable
#     _type_prefix = rhs._type_prefix;
            self._type_prefix = rhs._type_prefix
#     _scaling_coefficient = rhs._scaling_coefficient;
            self._scaling_coefficient = rhs._scaling_coefficient
#     for (unsigned int i=0; i<rhs._predicate_group_array.size(); i++){
            for ppp in rhs._predicate_group_array:
#         _predicate_group_array.push_back(new predicate_group_m_t(*(rhs._predicate_group_array[i])));
                self._predicate_group_array.append(PredicateGroupMT(rhs=ppp))
#     }
# }

# resource_m_t::~resource_m_t (){
#     for (unsigned int i=0; i<_predicate_group_array.size(); i++){
#         delete _predicate_group_array[i];
#     }
# }

# void resource_m_t::generate (const namespace_map & n_map, map<string, unsigned int> & id_cursor_map){
    def generate(self, n_map, id_cursor_map):
        """

        Args:

        Returns:

        """
#     if (id_cursor_map.find(_type_prefix)==id_cursor_map.end()){
        try:
            xxx = id_cursor_map[self._type_prefix]
        except KeyError:
#         id_cursor_map[_type_prefix] = 0;
            id_cursor_map[self._type_prefix] = 0
#     }
#     for (unsigned int id=id_cursor_map[_type_prefix]; id<(id_cursor_map[_type_prefix] + _scaling_coefficient); id++){
        for i in range(id_cursor_map[self._type_prefix], id_cursor_map[self._type_prefix] + self._scaling_coefficient):
#         string subject = "";
            subject = ''
#         subject.append("<");
            subject += '<'
#         subject.append(n_map.replace(_type_prefix));
            subject += n_map.replace(self._type_prefix)
#         subject.append(boost::lexical_cast<string>(id));
            subject += str(id)
#         subject.append(">");
            subject += '>'

#         for (vector<predicate_group_m_t*>::const_iterator itr2=_predicate_group_array.begin(); itr2!=_predicate_group_array.end(); itr2++){
            for itr2 in self._predicate_group_array:
#             predicate_group_m_t * predicate_group = *itr2;
                predicate_group = itr2
#             if (!predicate_group->_post_process){
                if not predicate_group._post_process:
#                 float draw = ((float) rand())/((float)RAND_MAX);
                    draw = random.uniform(0, 1)
#                 if (draw<=predicate_group->_gen_probability){
                    if draw <= predicate_group._gen_probability:
#                     for (vector<predicate_m_t*>::const_iterator itr3=predicate_group->_predicate_array.begin(); itr3!=predicate_group->_predicate_array.end(); itr3++){
                        for itr3 in predicate_group._predicate_array:
#                         predicate_m_t * predicate = *itr3;
                            predicate = itr3
#                         string triple_str = "";
                            triple_str = ''
#                         triple_str.append(subject);
                            triple_str += subject
#                         triple_str.append("\t");
                            triple_str += '\t'
#                         triple_str.append(predicate->generate(n_map));
                            triple_str += predicate.generate(n_map)

#                         int tab1_index = triple_str.find("\t");
                            tab1_index = triple_str.find('\t')
#                         int tab2_index = triple_str.find("\t", tab1_index+1);
                            tab2_index = triple_str.find('\t', tab1_index+1)

#                         //triple_lines.push_back(triple_st(triple_str.substr(0, tab1_index), triple_str.substr((tab1_index+1), (tab2_index-tab1_index-1)), triple_str.substr(tab2_index+1)));
#                         triple_st line (triple_str.substr(0, tab1_index), triple_str.substr((tab1_index+1), (tab2_index-tab1_index-1)), triple_str.substr(tab2_index+1));
                            line = TripleST(subject=triple_str[:tab1_index], predicate=triple_str[tab1_index+1:tab2_index], object0=triple_str[tab2_index + 1:])
#                         cout<<line<<" .\n";
                            print(line, ' .')
#                     }
#                 }
#             }
#         }
#     }
#     id_cursor_map[_type_prefix] += _scaling_coefficient;
        id_cursor_map[self._type_prefix] += self._scaling_coefficient
# }
        pass  # end of generate

# void resource_m_t::process_type_restrictions (const namespace_map & n_map, const type_map & t_map, const map<string, unsigned int> & id_cursor_map){
    def process_type_restrictions(self, n_map, t_map, id_cursor_map):
        """

        Args:

        Returns:

        """
#     unsigned int max_count = (id_cursor_map.find(_type_prefix))->second;
        max_count: int = id_cursor_map[self._type_prefix]
#     for (unsigned int id=0; id<max_count; id++){
        for i in range(max_count):
#         string subject = "";
            subject = ''
#         subject.append(n_map.replace(_type_prefix));
            subject += n_map.replace(self._type_prefix)
#         subject.append(boost::lexical_cast<string>(id));
            subject += str(id)

#         for (vector<predicate_group_m_t*>::const_iterator itr2=_predicate_group_array.begin(); itr2!=_predicate_group_array.end(); itr2++){
            for itr2 in self._predicate_group_array:
#             predicate_group_m_t * predicate_group = *itr2;
                predicate_group = itr2
#             if (predicate_group->_post_process && t_map.instanceof(subject, n_map.replace(*(predicate_group->_type_restriction)))){
                if predicate_group._post_process and t_map.instance_of(subject, n_map.replace(predicate_group._type_restriction)):
#                 float draw = ((float) rand())/((float)RAND_MAX);
                    draw = random.uniform(0, 1)
#                 if (draw<=predicate_group->_gen_probability){
                    if draw <= predicate_group._gen_probability:
#                     for (vector<predicate_m_t*>::const_iterator itr3=predicate_group->_predicate_array.begin(); itr3!=predicate_group->_predicate_array.end(); itr3++){
                        for itr3 in predicate_group._predicate_array:
#                         predicate_m_t * predicate = *itr3;
                            predicate = itr3
#                         string triple_str = "";
                            triple_str = ''
#                         triple_str.append("<");
                            triple_str += '<'
#                         triple_str.append(subject);
                            triple_str += subject
#                         triple_str.append(">");
                            triple_str += '>'
#                         triple_str.append("\t");
                            triple_str += '\t'
#                         triple_str.append(predicate->generate(n_map));
                            triple_str += predicate.generate(n_map)

#                         int tab1_index = triple_str.find("\t");
                            tab1_index: int = triple_str.find('\t')
#                         int tab2_index = triple_str.find("\t", tab1_index+1);
                            tab2_index: int = triple_str.find('\t', tab1_index+1)

#                         //triple_lines.push_back(triple_st(triple_str.substr(0, tab1_index), triple_str.substr((tab1_index+1), (tab2_index-tab1_index-1)), triple_str.substr(tab2_index+1)));
#                         triple_st line (triple_str.substr(0, tab1_index), triple_str.substr((tab1_index+1), (tab2_index-tab1_index-1)), triple_str.substr(tab2_index+1));
                            line = TripleST(subject=triple_str[:tab1_index], predicate=triple_str[tab1_index+1: tab2_index], object0=triple_str[tab2_index + 1:])
#                         cout<<line<<" .\n";
                            print(line, ' .')
#                     }
#                 }
#             }
#         }
#     }
# }
        pass  # end of

    @staticmethod
# resource_m_t * resource_m_t::parse (const string & line){
    def parse(line):
        """

        Args:

        Returns:

        """
#     bool scalable = true;
        scalable = True
#     string type_prefix ("not_defined");
        type_prefix = 'not_defined'
#     unsigned int scaling_coefficient = 0;
        scaling_coefficient: int = 0

#     stringstream parser(line);
        parser = parse_line(line)
#     int index = 0;
        index = 0
#     string token;
#     while (parser>>token){
        for token in parser:
#         switch (index){
#             case 0:{
            if index == 0:
#                 if (token.compare("<type>")!=0 && token.compare("<type*>")!=0 ){
                if token != '<type>' and token != '<type*>':
#                     cerr<<"[resource_m_t::parse()]\tExpecting <type> or <type*>..."<<"\n";
                    print('[resource_m_t::parse()]\tExpecting <type> or <type*>...')
#                     exit(0);
                    sys.exit(0)
#                 }
#                 if (token.compare("<type*>")==0){
                if token == '<type*>':
#                     scalable = false;
                    scalable = False
#                 }
#                 break;
#             }
#             case 1:{
            if index == 1:
#                 type_prefix = token;
                type_prefix = token
#                 break;
#             }
#             case 2:{
            if index == 2:
#                 scaling_coefficient = boost::lexical_cast<unsigned int>(token);
                scaling_coefficient = int(token)
#                 break;
#             }
#         }
#         index++;
            index += 1
#     }

#     if (index==3){
        if index == 3:
#         return new resource_m_t(scalable, type_prefix, scaling_coefficient);
            return ResourceMT(scalable=scalable, type_prefix=type_prefix, scaling_coefficient=scaling_coefficient)
#     } else {
        else:
#         cerr<<"[resource_m_t::parse()]\tExpecting 2 arguments..."<<"\n";
            print('[resource_m_t::parse()]\tExpecting 2 arguments...')
#         exit(0);
            sys.exit(0)
#     }

# }
        pass  # enf of parse
    pass  # end of class ResourceMT


# class namespace_map {
class NamespaceMap:
    """

    Attributes:

    """
#     public:
#         namespace_map();
#         ~namespace_map();

#         void insert (const namespace_m_t & namespace_declaration);
#         void insert (const string & alias, const string & prefix);
#         string lookup (const string & alias) const;
#         string replace (const string & content) const;

#         void to_str (vector<string> & lines) const;
# };

# namespace_map::namespace_map(){
    def __init__(self):
        """

        Args:

        Returns:

        """
#      private:
#         map<string, string> _index; // @suppress("Invalid template argument")
        self._index = {}
# }
        pass  # end of __init__

# namespace_map::~namespace_map(){

# }

# void namespace_map::insert (const namespace_m_t & namespace_declaration){
    def insert(self, namespace_declaration=None, alias=None, prefix=None):
        """

        Args:

        Returns:

        """
        if namespace_declaration is not None:
    #     if (_index.find(namespace_declaration._alias)==_index.end()){
            try:
                xxx = self._index[namespace_declaration._alias]
    #         _index.insert(pair<string, string>(namespace_declaration._alias, namespace_declaration._prefix));
    #     } else {
    #         cerr<<"[namespace_map::insert()] Warning: trying to insert an already existing namespace declaration..."<<"\n";
                print('[namespace_map::insert()] Warning: trying to insert an already existing namespace declaration...')
            except KeyError:
                self._index[namespace_declaration._alias] = namespace_declaration._prefix
    #     }
    # }
        pass  # end of insert

# void namespace_map::insert (const string & alias, const string & prefix){
        if alias is not None and prefix is not None:
#     if (_index.find(alias)==_index.end()){
            try:
                xxx = self._index[alias]
#         _index.insert(pair<string, string>(alias, prefix));
#     } else {
#         cerr<<"[namespace_map::insert()] Warning: trying to insert an already existing namespace declaration..."<<"\n";
                print('[namespace_map::insert()] Warning: trying to insert an already existing namespace declaration...')
            except KeyError:
                self._index[alias] = prefix
#     }
# }
        pass  # end of namespace_map::insert (const string & alias, const string & prefix)

# string namespace_map::lookup (const string & alias) const{
    def lookup(self, alias):
        """

        Args:

        Returns:

        """
#     map<string, string>::const_iterator f_it = _index.find(alias);
        try:
            f_it = self._index[alias]
#     if (f_it!=_index.end()){
#         return f_it->second;
            return f_it
#     } else {
        except KeyError:
#         cerr<<"[namespace_map::lookup()] Error: alias does not exist..."<<"\n";
            print('[namespace_map::lookup()] Error: alias does not exist...')
#         exit(0);
            sys.exit(0)
#     }
# }
        pass  # end of lookup

# string namespace_map::replace (const string & content) const{
    def replace(self, content):
        """

        Args:

        Returns:

        """
#     unsigned int pos = content.find_first_of(':');
        pos: int = content.find(':')
#     if (pos!=string::npos){
        if pos >= 0:
#         string alias = content.substr(0, pos);
            alias = content[:pos]
#         string suffix = content.substr(pos+1);
            suffix = content[pos+1:]
#         string result = lookup(alias);
            result = self.lookup(alias)
#         result.append(suffix);
            result += suffix
#         return result;
            return result
#     } else {
        else:
#         return content;
            return content
#     }
# }
        pass  # end of replace
    pass  # end of NamespaceMap


# struct association_m_t {
class AssociationMT:
    """

    Attributes:

    """

#     void init (string subject_type, string predicate, string object_type);

#     association_m_t (string subject_type, string predicate, string object_type);
#     association_m_t (string subject_type, string predicate, string object_type, unsigned int left_cardinality, unsigned int right_cardinality);
#     association_m_t (string subject_type, string predicate, string object_type, unsigned int left_cardinality, unsigned int right_cardinality, float left_cover);
#     association_m_t (string subject_type, string predicate, string object_type, unsigned int left_cardinality, unsigned int right_cardinality, float left_cover, DISTRIBUTION_TYPES::enum_t right_distribution);
#     association_m_t (string subject_type, string predicate, string object_type, unsigned int left_cardinality, unsigned int right_cardinality, float left_cover, DISTRIBUTION_TYPES::enum_t right_distribution, const string * subject_type_restriction, const string * object_type_restriction);
#     ~association_m_t ();

#     void generate (const namespace_map & n_map, type_map & t_map, const map<string, unsigned int> & id_cursor_map); // @suppress("Invalid template argument")
#     void process_type_restrictions (const namespace_map & n_map, const type_map & t_map, const map<string, unsigned int> & id_cursor_map); // @suppress("Invalid template argument")

#     static association_m_t * parse (const map<string, unsigned int> & id_cursor_map, const string & line); // @suppress("Invalid template argument")
# };

# void association_m_t::init (string subject_type, string predicate, string object_type){
    def init(self, subject_type: str, predicate: str, object_type: str):
        """

        Args:

        Returns:

        """
#     _post_process = false;
        self._post_process = False
#     _subject_type = subject_type;
        self._subject_type = subject_type
#     _predicate = predicate;
        self._predicate = predicate
#     _object_type = object_type;
        self._object_type = object_type
#     _left_cardinality = 1;
        self._left_cardinality = 1
#     _right_cardinality = 1;
        self._right_cardinality = 1
#     _right_cardinality_distribution = DISTRIBUTION_TYPES::UNDEFINED;
        self._right_cardinality_distribution = DistributionTypes.UNDEFINED
#     _left_cover = 1.0;
        self._left_cover = 1.0
#     _right_distribution = DISTRIBUTION_TYPES::UNIFORM;
        self._right_distribution = DistributionTypes.UNIFORM
#     _subject_type_restriction = NULL;
        self._subject_type_restriction = None
#     _object_type_restriction = NULL;
        self._object_type_restriction = None
# }
        pass  # end of init

# association_m_t::association_m_t (string subject_type, string predicate, string object_type){
    def __init__(self, subject_type, predicate, object_type, left_cardinality=None, right_cardinality=None, left_cover=None,
                 right_distribution=None, subject_type_restriction=None, object_type_restriction=None):
#     bool                            _post_process;
        self._post_process = False
#     string                          _subject_type;
        self._subject_type = ''
#     string                          _predicate;
        self._predicate = ''
#     string                          _object_type;
        self._object_type = ''

#     string *                        _subject_type_restriction;
        self._subject_type_restriction = ''
#     string *                        _object_type_restriction;
        self._object_type_restriction = ''

#     unsigned int                    _left_cardinality;
        self._left_cardinality = 0
#     unsigned int                    _right_cardinality;
        self._right_cardinality = 0
#     DISTRIBUTION_TYPES::enum_t      _right_cardinality_distribution;
        self._right_cardinality_distribution = None

#     float                           _left_cover;
        self._left_cover = 0.0
#     DISTRIBUTION_TYPES::enum_t      _right_distribution;
        self._right_distribution = None
#     init (subject_type, predicate, object_type);
        self.init(subject_type, predicate, object_type)
# }

# association_m_t::association_m_t (string subject_type, string predicate, string object_type, unsigned int left_cardinality, unsigned int right_cardinality){
        if left_cardinality is not None and right_cardinality is not None:
#     init (subject_type, predicate, object_type);
#     _left_cardinality = left_cardinality;
            self._left_cardinality = left_cardinality
#     _right_cardinality = right_cardinality;
            self._right_cardinality = right_cardinality
# }

# association_m_t::association_m_t (string subject_type, string predicate, string object_type, unsigned int left_cardinality, unsigned int right_cardinality,
# float left_cover){
        if left_cover is not None:
#     init (subject_type, predicate, object_type);
#     _left_cardinality = left_cardinality;
#     _right_cardinality = right_cardinality;
#     _left_cover = left_cover;
            self._left_cover = left_cover
# }

# association_m_t::association_m_t (string subject_type, string predicate, string object_type, unsigned int left_cardinality, unsigned int right_cardinality,
# float left_cover, DISTRIBUTION_TYPES::enum_t right_distribution){
        if right_distribution is not None:
#     init (subject_type, predicate, object_type);
#     _left_cardinality = left_cardinality;
#     _right_cardinality = right_cardinality;
#     _left_cover = left_cover;
#     _right_distribution = right_distribution;
            self._right_distribution = right_distribution
# }

# association_m_t::association_m_t (
#     string subject_type,
#     string predicate,
#     string object_type,
#     unsigned int left_cardinality,
#     unsigned int right_cardinality,
#     float left_cover,
#     DISTRIBUTION_TYPES::enum_t right_distribution,
#     const string * subject_type_restriction,
#     const string * object_type_restriction)
        if subject_type_restriction is not None and object_type_restriction is not None:
#     {
#         init (subject_type, predicate, object_type);
#         _post_process = true;
            self._post_process = True
#         _left_cardinality = left_cardinality;
#         _right_cardinality = right_cardinality;
#         _left_cover = left_cover;
#         _right_distribution = right_distribution;
#         if (subject_type_restriction!=NULL){
#             _subject_type_restriction = new string(*subject_type_restriction);
            self._subject_type_restriction = str(subject_type_restriction)
#         }
#         if (object_type_restriction!=NULL){
#             _object_type_restriction = new string(*object_type_restriction);
            self._object_type_restriction = str(object_type_restriction)
#         }
# }
        pass  # end of

        """
    
        Args:
    
        Returns:
    
        """

    # association_m_t::~association_m_t (){
#     delete _subject_type_restriction;
#     delete _object_type_restriction;
# }

# void association_m_t::generate (const namespace_map & n_map, type_map & t_map, const map<string, unsigned int> & id_cursor_map){
    def generate(self, n_map, t_map: TypeMap, id_cursor_map):
        """

        Args:

        Returns:

        """
#     if (id_cursor_map.find(_subject_type)==id_cursor_map.end()){
        try:
            xxx = id_cursor_map[self._subject_type]
        except KeyError:
#         cerr<<"[association_m_t::parse()] Error: association cannot be defined over undefined resource '"<<_subject_type<<"'..."<<"\n";
            print(f'[association_m_t::parse()] Error: association cannot be defined over undefined resource "{self._subject_type}"...')
#         exit(0);
            sys.exit(0)
#     }
#     if (id_cursor_map.find(_object_type)==id_cursor_map.end()){
        try:
            xxx = id_cursor_map[self._object_type]
        except KeyError:
#         cerr<<"[association_m_t::parse()] Error: association cannot be defined over undefined resource '"<<_object_type<<"'..."<<"\n";
            print(f'[association_m_t::parse()] Error: association cannot be defined over undefined resource "{self._object_type}"...')
#         exit(0);
            sys.exit(0)
#     }

#     model::clear_zipfian_cache();
        Model.clear_zipfian_cache()

#     if (!_post_process){
        if not self._post_process:
#         unsigned int left_instance_count = id_cursor_map.find(_subject_type)->second;
            left_instance_count = id_cursor_map[self._subject_type]
#         unsigned int right_instance_count = id_cursor_map.find(_object_type)->second;
            right_instance_count = id_cursor_map[self._object_type]
#         unordered_set<unsigned int> mapped_instances;
            mapped_instances = set()
#
#         boost::posix_time::ptime t1 (bpt::microsec_clock::universal_time());

#         for (unsigned int left_id=0; left_id<left_instance_count; left_id++){
            for left_id in range(left_instance_count):
#             float pr = ((float) rand()) / ((float) RAND_MAX);
                pr = random.uniform(0, 1)
#             if (pr<=_left_cover){
                if pr <= self._left_cover:
#                 unsigned int right_size = _right_cardinality;
                    right_size = self._right_cardinality
#                 if (_right_cardinality_distribution!=DISTRIBUTION_TYPES::UNDEFINED){
                    if self._right_cardinality_distribution != DistributionTypes.UNDEFINED:
#                     right_size = round((double) right_size * model::generate_random(_right_cardinality_distribution));
                        right_size = round(right_size * Model.generate_random(self._right_cardinality_distribution))
#                     right_size = (right_size > _right_cardinality) ? _right_cardinality : right_size;
                        if right_size > self._right_cardinality:
                            right_size = self._right_cardinality
#                 }
#                 for (unsigned int j=0; j<right_size; j++){
                    for j in range(right_size):
#                     unsigned int loop_counter = 0;
                        loop_counter = 0
#                     unsigned int right_id = 0;
                        right_id = 0
#                     do {
                        while True:
#                         double r_value = model::generate_random(_right_distribution, right_instance_count);
                            r_value = Model.generate_random(self._right_distribution, right_instance_count)
#                         right_id = round(r_value * right_instance_count);
                            right_id = round(r_value * right_instance_count)
#                         right_id = (right_id>=right_instance_count) ? (right_instance_count-1) : right_id;
                            if right_id >= right_instance_count:
                                right_id = right_instance_count - 1
#                         loop_counter++;
                            loop_counter += 1
#                     } while (mapped_instances.find(right_id)!=mapped_instances.end() && loop_counter<MAX_LOOP_COUNTER);
                            if right_id not in mapped_instances:
                                break
                            if loop_counter >= MAX_LOOP_COUNTER:
                                break
#                     if (loop_counter<MAX_LOOP_COUNTER){
                        if loop_counter < MAX_LOOP_COUNTER:
#                         if (_left_cardinality==1){
                            if self._left_cardinality == 1:
#                             mapped_instances.insert(right_id);
                                mapped_instances.add(right_id)
#                         }
#                         string subject(""), predicate (""), object(""), triple ("");
                            subject = ''
                            predicate = ''
                            object0 = ''
                            triple = ''

#                         // FIXME:: You need to add replace-command...
#                         subject.append(n_map.replace(_subject_type));
                            subject += n_map.replace(self._subject_type)
#                         subject.append(boost::lexical_cast<string>(left_id));
                            subject += str(left_id)

#                         object.append(n_map.replace(_object_type));
                            object0 += n_map.replace(self._object_type)
#                         object.append(boost::lexical_cast<string>(right_id));
                            object0 += str(right_id)

#                         predicate.append(n_map.replace(_predicate));
                            predicate += n_map.replace(self._predicate)

#                         string subject_str(""), predicate_str (""), object_str("");
                            subject_str = ''
                            predicate_str = ''
                            object_str = ''
#                         subject_str.append("<");
                            subject_str += '<'
#                         subject_str.append(subject);
                            subject_str += subject
#                         subject_str.append(">");
                            subject_str += '>'

#                         predicate_str.append("<");
                            predicate_str += '<'
#                         predicate_str.append(predicate);
                            predicate_str += predicate
#                         predicate_str.append(">");
                            predicate_str += '>'

#                         object_str.append("<");
                            object_str += '<'
#                         object_str.append(object);
                            object_str += object0
#                         object_str.append(">");
                            object_str += '>'

#                         //triple_lines.push_back(triple_st(subject_str, predicate_str, object_str));
#                         triple_st line (subject_str, predicate_str, object_str);
                            line = TripleST(subject=subject_str, predicate=predicate_str, object0=object_str)
#                         cout<<line<<" .\n";
                            print(line, ' .')

#                         // Save type assertions...
#                         if (predicate.compare("http://www.w3.org/1999/02/22-rdf-syntax-ns#type")==0){
                            if predicate == 'http://www.w3.org/1999/02/22-rdf-syntax-ns#type':
#                             t_map.insert(subject, object);
                                t_map.insert(subject, object0)
#                         }
#                     } else {
                        else:
#                         //cout<<"[association_m_t::generate] Warning:: failed to greedily satisfy cardinality constraints..."<<"\n";
#                         //cout<<"[association_m_t::generate] Ignoring association "
#                             //<<_subject_type<<left_id<<"-->"
#                             //<<_object_type<<right_id<<" after"<<MAX_LOOP_COUNTER<<" trials..."<<"\n";
                            pass
#                     }
#                 }
#             }
#         }

#         boost::posix_time::ptime t2 (bpt::microsec_clock::universal_time());
#         //cerr    << "[association-generation]" << " " << (t2-t1).total_microseconds() << " "
#         //        << _subject_type << " " << _predicate << " " << _object_type << " "
#         //        << _left_cardinality << " " << _right_cardinality << " "
#         //        << "\n";
#     }
# }
        pass  # end of

# void association_m_t::process_type_restrictions (const namespace_map & n_map, const type_map & t_map, const map<string, unsigned int> & id_cursor_map){
    def process_type_restrictions(self, n_map, t_map, id_cursor_map):
        """

        Args:

        Returns:

        """
        right_instance_count = 0
#     if (id_cursor_map.find(_subject_type)==id_cursor_map.end()){
        try:
            xxx = id_cursor_map[self._subject_type]
        except KeyError:
#         cerr<<"[association_m_t::parse()] Error: association cannot be defined over undefined resource '"<<_subject_type<<"'..."<<"\n";
            print(f'[association_m_t::parse()] Error: association cannot be defined over undefined resource "{self._subject_type}"...')
#         exit(0);
            sys.exit(0)
#     }
#     if (id_cursor_map.find(_object_type)==id_cursor_map.end()){
        try:
            xxx = id_cursor_map[self._object_type]
        except KeyError:
#         cerr<<"[association_m_t::parse()] Error: association cannot be defined over undefined resource '"<<_object_type<<"'..."<<"\n";
            print(f'[association_m_t::parse()] Error: association cannot be defined over undefined resource "{self._object_type}"...')
#         exit(0);
            sys.exit(0)
#     }

#     if (_post_process){
        if self._post_process:
#         unsigned int left_instance_count = id_cursor_map.find(_subject_type)->second;
            left_instance_count = id_cursor_map[self._subject_type]
#         vector<string> * restricted_right_instances = NULL;
            restricted_right_instances = []
#         if (_object_type_restriction!=NULL){
            if self._object_type_restriction:
#             restricted_right_instances = t_map.get_instances(n_map.replace(_object_type), n_map.replace(*_object_type_restriction));
                restricted_right_instances = t_map.get_instances(n_map.replace(self._object_type), n_map.replace(self._object_type_restriction))
#         } else {
            else:
#             restricted_right_instances = new vector<string>();
                restricted_right_instances = []
#             unsigned int right_instance_count = id_cursor_map.find(_object_type)->second;
                right_instance_count = id_cursor_map[self._object_type]
#             for (unsigned int right_id=0; right_id<right_instance_count; right_id++){
                for right_id in range(right_instance_count):
#                 string instance = n_map.replace(_object_type);
                    instance = n_map.replace(self._object_type)
#                 instance.append(boost::lexical_cast<string>(right_id));
                    instance += str(right_id)
#                 restricted_right_instances->push_back(instance);
                    restricted_right_instances.append(instance)
#             }
#         }
#         if (restricted_right_instances!=NULL){
            if restricted_right_instances:
#             unsigned int right_instance_count = restricted_right_instances->size();
                right_instance_coiunt = len(restricted_right_instances)
#             set<string> mapped_instances;
                mapped_instances = set()
#             for (unsigned int left_id=0; left_id<left_instance_count; left_id++){
                for left_id in range(left_instance_count):
#                 string subject="";
                    subject = ''
#                 subject.append(n_map.replace(_subject_type));
                    subject += n_map.replace(self._subject_type)
#                 subject.append(boost::lexical_cast<string>(left_id));
                    subject += str(left_id)
#                 if (_subject_type_restriction==NULL || t_map.instanceof(subject, n_map.replace(*_subject_type_restriction))){
                    if self._subject_type_restriction == None or t_map.instance_of(subject, n_map.replace(self._subject_type_restriction)):
#                     float pr = ((float) rand()) / ((float) RAND_MAX);
                        pr = random.uniform(0, 1)
#                     if (pr<=_left_cover){
                        if pr <= self._left_cover:
#                         unsigned int right_size = _right_cardinality;
                            right_size = self._right_cardinality
#                         if (_right_cardinality_distribution!=DISTRIBUTION_TYPES::UNDEFINED){
                            if self._right_cardinality_distribution != DistributionTypes.UNDEFINED:
#                             right_size = round((double) right_size * model::generate_random(_right_cardinality_distribution));
                                right_size = round(right_size * Model.generate_random(self._right_cardinality_distribution))
#                             right_size = (right_size > _right_cardinality) ? _right_cardinality : right_size;
                                if right_size > self._right_cardinality:
                                    right_size = self._right_cardinality
#                         }
#                         for (unsigned int j=0; j<right_size; j++){
                            for j in range(right_size):
#                             string predicate="", object="", triple="";
                                predicate = ''
                                object0 = ''
                                triple = ''
#                             unsigned int loop_counter = 0;
                                loop_counter = 0
#                             do {
                                while True:
#                                 double r_value = model::generate_random(_right_distribution, right_instance_count);
                                    r_value = Model.generate_random(self._right_distribution, right_instance_count)
#                                 unsigned int right_index = round(r_value * right_instance_count);
                                    right_index = round(r_value * right_instance_count)
#                                 right_index = (right_index>=right_instance_count) ? (right_instance_count-1) : right_index;
                                    if right_index >= right_instance_count:
                                        right_index = right_instance_count - 1
#                                 object = (*restricted_right_instances)[right_index];
                                    object0 = restricted_right_instances[right_index]
#                                 loop_counter++;
                                    loop_counter += 1
#                             } while (mapped_instances.find(object)!=mapped_instances.end() && loop_counter<MAX_LOOP_COUNTER);
                                    if object0 in mapped_instances:
                                        break
                                    if loop_counter >= MAX_LOOP_COUNTER:
                                        break
#                             if (loop_counter<MAX_LOOP_COUNTER){
                                if loop_counter < MAX_LOOP_COUNTER:
#                                 if (_left_cardinality==1){
                                    if self._left_cardinality == 1:
#                                     mapped_instances.insert(object);
                                        mapped_instances.add(object0)
#                                 }

#                                 predicate.append(n_map.replace(_predicate));
                                    predicate += n_map.replaces(self._predicate)

#                                 string subject_str = "", predicate_str = "", object_str = "";
                                    subject_str = ''
                                    predicate_str = ''
                                    object_str = ''
#                                 subject_str.append("<");
                                    subject_str += '<'
#                                 subject_str.append(subject);
                                    subject_str += subject
#                                 subject_str.append(">");
                                    subject_str += '>'
#                                 predicate_str.append("<");
                                    predicate_str += '<'
#                                 predicate_str.append(predicate);
                                    predicate_str += predicate
#                                 predicate_str.append(">");
                                    predicate_str += '>'
#                                 object_str.append("<");
                                    object_str += '<'
#                                 object_str.append(object);
                                    object_str += object0
#                                 object_str.append(">");
                                    object_str += '>'

#                                 //triple_lines.push_back(triple_st(subject_str, predicate_str, object_str));
#                                 triple_st line (subject_str, predicate_str, object_str);
                                    line = TripleST(subject=subject_str, predicate=predicate_str, object0=object_str)
#                                 cout<<line<<" .\n";
                                    print(line, ' .')
#                             }
#                         }
#                     }
#                 }
#             }
#             delete restricted_right_instances;
#         }
#     }
# }
        pass  # end of

    @staticmethod
# association_m_t * association_m_t::parse (const map<string, unsigned int> & id_cursor_map, const string & line){
    def parse(line):
        """

        Args:

        Returns:

        """
#     string subject_type ("");
        subject_type = ''
#     string predicate ("");
        predicate = ''
#     string object_type ("");
        object_type = ''
#     unsigned left_cardinality = 1;
        left_cardinality = 1
#     unsigned right_cardinality = 1;
        right_cardinality = 1
#     DISTRIBUTION_TYPES::enum_t right_cardinality_distribution = DISTRIBUTION_TYPES::UNDEFINED;
        right_cardinality_distribution = DistributionTypes.UNDEFINED
#     float left_cover = 1.0;
        left_cover = 1.0
#     DISTRIBUTION_TYPES::enum_t right_distribution = DISTRIBUTION_TYPES::UNIFORM;
        right_distribution = DistributionTypes.UNIFORM
#     string * subject_type_restriction = NULL;
        subject_type_restriction = None
#     string * object_type_restriction = NULL;
        object_type_restriction = None

#     stringstream parser(line);
        parser = parse_line(line)
#     int index = 0;
        index = 0
#     string token;
#     while (parser>>token){
        for token in parser:
#         switch (index){
#             case 0: {
            if index == 0:
#                 if (token.compare("#association")!=0){
                if token != '#association':
#                     cerr<<"[association_m_t::parse()]\tExpecting #association..."<<"\n";
                    print('[association_m_t::parse()]\tExpecting #association...')
#                     exit(0);
                    sys.exit(0)
#                 }
#                 break;
#             }
#             case 1: {
            elif index == 1:
#                 subject_type = token;
                subject_type = token
#                 break;
#             }
#             case 2: {
            elif index == 2:
#                 predicate = token;
                predicate = token
#                 break;
#             }
#             case 3: {
            elif index == 3:
#                 object_type = token;
                object_type = token
#                 break;
#             }
#             case 4: {
            elif index == 4:
#                 left_cardinality = boost::lexical_cast<unsigned int>(token);
                left_cardinality = int(token)
#                 break;
#             }
#             case 5: {
            elif index == 5:
#                 if (token.find("[uniform]")!=string::npos || token.find("[UNIFORM]")!=string::npos){
                if token.find('[uniform]') >= 0 or token.find('[UNIFORM]') >= 0:
#                     right_cardinality_distribution = DISTRIBUTION_TYPES::UNIFORM;
                    right_cardinality_distribution = DistributionTypes.UNIFORM
#                     token = token.substr(0, token.find_first_of('['));
                    token = token[:token.find('[')]
#                 } else if (token.find("[normal]")!=string::npos || token.find("[NORMAL]")!=string::npos){
                elif token.find('[normal]') >= 0 or token.find('[NORMAL]') >= 0:
    #                     right_cardinality_distribution = DISTRIBUTION_TYPES::NORMAL;
                    right_cardinality_distribution = DistributionTypes.NORMAL
#                     token = token.substr(0, token.find_first_of('['));
                    token = token[:token.find('[')]
#                 }
#                 right_cardinality = boost::lexical_cast<unsigned int>(token);
                right_cardinality = int(token)
#                 break;
#             }
#             case 6: {
            elif index == 6:
#                 left_cover = boost::lexical_cast<float>(token);
                left_cover = float(token)
#                 break;
#             }
#             case 7: {
            elif index == 7:
#                 if (token.compare("uniform")==0 || token.compare("UNIFORM")==0){
                if token.find('[uniform]') >= 0 or token.find('[UNIFORM]') >= 0:
#                     right_distribution = DISTRIBUTION_TYPES::UNIFORM;
                    right_distribution = DistributionTypes.UNIFORM
#                 } else if (token.compare("normal")==0 || token.compare("NORMAL")==0){
                elif token.find('[normal]') >= 0 or token.find('[NORMAL]') >= 0:
#                     right_distribution = DISTRIBUTION_TYPES::NORMAL;
                    right_distribution = DistributionTypes.NORMAL
#                 } else if (token.compare("zipfian")==0 || token.compare("ZIPFIAN")==0){
                elif token.find('[normal]') >= 0 or token.find('[NORMAL]') >= 0:
#                     right_distribution = DISTRIBUTION_TYPES::ZIPFIAN;
                    right_distribution = DistributionTypes.ZIPFIAN
#                 }
#                 break;
#             }
#             case 8: {
            elif index == 8:
#                 if (!boost::starts_with(token, "@")){
                if not token.startswith('@'):
#                     cerr<<"[association_m_t::parse()]\tExpecting '@' qualifier before type restriction..."<<"\n";
                    print('[association_m_t::parse()]\tExpecting '@' qualifier before type restriction...')
#                     exit(0);
                    sys.exit(0)
#                 }
#                 if (token.compare("@null")!=0 && token.compare("@NULL")!=0){
                if token.find('@null') < 0 and token.find('@NULL') < 0:
#                     subject_type_restriction = new string(token.substr(1));
                    subject_type_restriction = str(token[1:])
#                 }
#                 break;
#             }
#             case 9: {
            elif index == 9:
#                 if (!boost::starts_with(token, "@")){
                if not token.startswith('@'):
#                     cerr<<"[association_m_t::parse()]\tExpecting '@' qualifier before type restriction..."<<"\n";
                    print('[association_m_t::parse()]\tExpecting '@' qualifier before type restriction...')
#                     exit(0);
                    sys.exit(0)
#                 }
#                 if (token.compare("@null")!=0 && token.compare("@NULL")!=0){
                if token.find('@null') < 0 and token.find('@NULL') < 0:
#                     object_type_restriction = new string(token.substr(1));
                    object_type_restriction = str(token[:1])
#                 }
#                 break;
#             }
#         }
#         index++;
            index += 1
#     }
#     association_m_t * result = NULL;
        result = None
#     if (index==4){
        if index == 4:
#         result = new association_m_t (subject_type, predicate, object_type);
            result = AssociationMT(subject_type, predicate, object_type)
#     } else if (index==6){
        elif index == 6:
#         result = new association_m_t (subject_type, predicate, object_type, left_cardinality, right_cardinality);
            result = AssociationMT(subject_type, predicate, object_type, left_cardinality=left_cardinality, right_cardinality=right_cardinality)
#     } else if (index==7){
        elif index == 7:
#         result = new association_m_t (subject_type, predicate, object_type, left_cardinality, right_cardinality, left_cover);
            result = AssociationMT(subject_type, predicate, object_type, left_cardinality=left_cardinality, right_cardinality=right_cardinality, left_cover=left_cover)
#     } else if (index==8){
        elif index == 8:
#         result = new association_m_t (subject_type, predicate, object_type, left_cardinality, right_cardinality, left_cover, right_distribution);
            result = AssociationMT(subject_type, predicate, object_type, left_cardinality=left_cardinality, right_cardinality=right_cardinality, left_cover=left_cover, right_distribution=right_distribution)
#     } else if (index==10){
        elif index == 10:
#         result = new association_m_t (subject_type, predicate, object_type, left_cardinality, right_cardinality, left_cover, right_distribution, subject_type_restriction, object_type_restriction);
            result = AssociationMT(subject_type, predicate, object_type, left_cardinality=left_cardinality, right_cardinality=right_cardinality, left_cover=left_cover, right_distribution=right_distribution, subject_type_restriction=subject_type_restriction, object_type_restriction=object_type_restriction)
        #     } else {
        else:
#         cerr<<"[association_m_t::parse()]\tExpecting 3, 5, 6, 7 or 9 arguments..."<<"\n";
            print('[association_m_t::parse()]\tExpecting 3, 5, 6, 7 or 9 arguments...')
#         exit(0);
            sys.exit(0)
#     }
#     result->_right_cardinality_distribution = right_cardinality_distribution;
        result.right_cardinality_distribution = right_cardinality_distribution
#     delete subject_type_restriction;
#     delete object_type_restriction;
#     return result;
        return result
# }
        pass # end of parse
    pass  # end of AssociationMT


# struct mapping_m_t {
class MappingMT:
    """

    Attributes:

    """
#     void init (const string & var_name, LITERAL_TYPES::enum_t literal_type);
#     void init (const string & var_name, const string & resource_type);

#     mapping_m_t (const string & var_name, LITERAL_TYPES::enum_t literal_type);
#     mapping_m_t (const string & var_name, LITERAL_TYPES::enum_t literal_type, DISTRIBUTION_TYPES::enum_t distribution_type);
#     mapping_m_t (const string & var_name, LITERAL_TYPES::enum_t literal_type, DISTRIBUTION_TYPES::enum_t distribution_type, const string & range_min, const string & range_max);
#     mapping_m_t (const string & var_name, const string & resource_type);
#     mapping_m_t (const string & var_name, const string & resource_type, const string & type_restriction);
#     mapping_m_t (const string & var_name, const string & resource_type, DISTRIBUTION_TYPES::enum_t distribution_type);
#     mapping_m_t (const string & var_name, const string & resource_type, const string & type_restriction, DISTRIBUTION_TYPES::enum_t distribution_type);

#     mapping_m_t (const mapping_m_t & rhs);
#     ~mapping_m_t ();

#     string generate (const model & mdl, const query_template_m_t & q_template);
#     string generate (const model & mdl, const query_template_m_t & q_template, unsigned int & instance_count);

#     static mapping_m_t * parse (const string & line);
# };
# void mapping_m_t::init (const string & var_name, LITERAL_TYPES::enum_t literal_type){
    def init(self, var_name='', literal_type=None, resource_type=None):
        """

        Args:

        Returns:

        """
        if literal_type is not None:
    #     _var_name = var_name;
            self._var_name = var_name
    #     _is_literal_type = true;
            self._is_literal_type = True
    #     _literal_type = literal_type;
            self._literal_type = literal_type
    #     _resource_type = "";
            self._resource_type = ''
    #     _type_restriction = NULL;
            self._type_restriction = None
    #     _distribution_type = DISTRIBUTION_TYPES::UNIFORM;
            self._distribution_type = DistributionTypes.UNIFORM

    #     switch (_literal_type){
    #         case LITERAL_TYPES::INTEGER:{
            if self._literal_type == LiteralTypes.INTEGER:
    #             _range_min = boost::lexical_cast<string>(numeric_limits<unsigned short>::min());
                self._range_min = -1000000
    #             _range_max = boost::lexical_cast<string>(numeric_limits<unsigned short>::max());
                self._range_max = 1000000
    #             break;
    #         }
    #         case LITERAL_TYPES::STRING:
    #         case LITERAL_TYPES::NAME:{
            elif self._literal_type == LiteralTypes.STRING or self._literal_type == LiteralTypes.NAME:
    #             _range_min = string("A");
                self._range_min = 'A'
    #             _range_max = string("z");
                self._range_max = 'z'
    #             break;
    #         }
    #         case LITERAL_TYPES::DATE:{
            elif self._literal_type == LiteralTypes.DATE:
    #             boost::posix_time::ptime cur_time (boost::posix_time::second_clock::local_time());
                cur_time = datetime.time()
    #             boost::gregorian::date cur_date = cur_time.date();
                cur_date = datetime.date()
    #             _range_min = string("1970-01-01");
                self._range_min = '1970-01-01'
    #             _range_max = boost::gregorian::to_iso_extended_string(cur_date);
                self._range_max = str(cur_date)
    #             break;
    #         }
    #     }
    #     _dynamic_model_name = "";
            self._dynamic_model_name = ''
    # }

# void mapping_m_t::init (const string & var_name, const string & resource_type){
        if resource_type is not None:

#     _var_name = var_name;
            self._var_name = var_name
#     _is_literal_type =false;
            self._is_literal_type =False
#     _literal_type = LITERAL_TYPES::UNDEFINED;
            self._literal_type = LiteralTypes.UNDEFINED
#     _resource_type = resource_type;
            self._resource_type = resource_type
#     _type_restriction = NULL;
            self._type_restriction = None
#     _distribution_type = DISTRIBUTION_TYPES::UNIFORM;
            self._distribution_type = DistributionTypes.UNIFORM
#     _range_min = "";
            self._range_min = ''
#     _range_max = "";
            self._range_max = ''
#     _dynamic_model_name = "";
            self._dynamic_model_name = ''
# }

# mapping_m_t::mapping_m_t (const string & var_name, LITERAL_TYPES::enum_t literal_type){
    def __init__(self, var_name='', literal_type=None, distribution_type=None, range_min=None, range_max=None, resource_type=None, type_restriction=None, rhs=None):
        """

        Args:

        Returns:

        """
        #     string                      _var_name;
        #     bool                        _is_literal_type;
        #     LITERAL_TYPES::enum_t       _literal_type;
        #     string                      _resource_type;
        #     string *                     _type_restriction;
        #     DISTRIBUTION_TYPES::enum_t  _distribution_type;
        #     string                      _range_min;
        #     string                      _range_max;
        #     string                      _dynamic_model_name;
        self._var_name = ''
        self._is_literal_type = False
        self._literal_type = None
        self._resource_type = ''
        self._type_restriction = ''
        self._distribution_type = None
        self._range_min = ''
        self._range_max = ''
        self._dynamic_model_name = ''

        if literal_type is not None:
#     init (var_name, literal_type);
            self.init(var_name=var_name, literal_type=literal_type)
# }

# mapping_m_t::mapping_m_t (const string & var_name, LITERAL_TYPES::enum_t literal_type, DISTRIBUTION_TYPES::enum_t distribution_type){
#     init (var_name, literal_type);
        if distribution_type is not None:
#     _distribution_type = distribution_type;
            self._distribution_type = distribution_type
# }

# mapping_m_t::mapping_m_t (const string & var_name, LITERAL_TYPES::enum_t literal_type, DISTRIBUTION_TYPES::enum_t distribution_type, const string & range_min, const string & range_max){
        if range_min is not None and range_max is not None:
#     init (var_name, literal_type);
#     _distribution_type = distribution_type;
#     _range_min = range_min;
            self._range_min = range_min
#     _range_max = range_max;
            self._range_max = range_max
# }

# mapping_m_t::mapping_m_t (const string & var_name, const string & resource_type){
        if resource_type is not None:
#     init (var_name, resource_type);
            self.init(var_name=var_name, resource_type=resource_type)
# }

# mapping_m_t::mapping_m_t (const string & var_name, const string & resource_type, const string & type_restriction){
        if type_restriction is not None:
#     init (var_name, resource_type);
#     _type_restriction = new string (type_restriction);
            self._type_restriction = str(type_restriction)
# }

# mapping_m_t::mapping_m_t (const string & var_name, const string & resource_type, DISTRIBUTION_TYPES::enum_t distribution_type){
#     init (var_name, resource_type);
#     _distribution_type = distribution_type;
# }

# mapping_m_t::mapping_m_t (const string & var_name, const string & resource_type, const string & type_restriction, DISTRIBUTION_TYPES::enum_t distribution_type){
#     init (var_name, resource_type);
#     _type_restriction = new string (type_restriction);
#     _distribution_type = distribution_type;
# }

# mapping_m_t::mapping_m_t (const mapping_m_t & rhs){
        if rhs is not None:
#     _var_name = rhs._var_name;
            self._var_name = rhs._var_name
#     _is_literal_type = rhs._is_literal_type;
            self._is_literal_type = rhs._is_literal_type
#     _literal_type = rhs._literal_type;
            self._literal_type = rhs._literal_type
#     _resource_type = rhs._resource_type;
            self._resource_type = rhs._resource_type
#     _type_restriction = new string (*(rhs._type_restriction));
            self._type_restriction = str(rhs._type_restriction)
#     _distribution_type = rhs._distribution_type;
            self._distribution_type = rhs._distribution_type
#     _range_min = rhs._range_min;
            self._range_min = rhs._range_min
#     _range_max = rhs._range_max;
            self._range_max = rhs._range_max
#     _dynamic_model_name = rhs._dynamic_model_name;
            self._dynamic_model_name = rhs._dynamic_model_name
# }

# mapping_m_t::~mapping_m_t (){
#     delete _type_restriction;
# }

# string mapping_m_t::generate (const model & mdl, const query_template_m_t & q_template){
    def generate(self, mdl, q_template, instance_count=None):
        """

        Args:

        Returns:

        """
        if instance_count is None:
    #     unsigned int count = 0;
            count = 0
    #     return generate(mdl, q_template, count);
            return self.generate(mdl, q_template, instance_count=count)
    # }

# string mapping_m_t::generate (const model & mdl, const query_template_m_t & q_template, unsigned int & instance_count){
        else:
#     if (_is_literal_type){
            if self._is_literal_type:
#         string result = "";
                result = ''
#         result.append(model::generate_literal(_literal_type, _distribution_type, _range_min, _range_max));
                result += Model.generate_literal(self._literal_type, self._distribution_type, self._range_min, self._range_max)
#         return result;
                return result
#     } else {
            else:
#         if (_type_restriction==NULL){
                if self._type_restriction is None:
#             string result = "";
                    result = ''
#             instance_count = mdl._id_cursor_map.find(_resource_type)->second;
                    instance_count = mdl._id_cursor_map[self._resource_type]

#             unsigned int id = 0;
                    id0 = 0
#             if (_distribution_type==DISTRIBUTION_TYPES::DYNAMIC){
                    if self._distribution_type == DistributionTypes.DYNAMIC:
#                 map<string, pair<volatility_gen*, float> >::const_iterator itr = q_template._volatility_table.find(_dynamic_model_name);
                        try:
                            itr = q_template._volatility_table[self._dynamic_model_name]
#                 if (itr==q_template._volatility_table.end()){
                        except KeyError:
#                     cerr << "[mapping_m_t::generate()]\tDynamic model " << _dynamic_model_name << " does not exist..." << "\n";
                            print(f'[mapping_m_t::generate()]\tDynamic model {self._dynamic_model_name} does not exist...')
#                     exit(0);
                            sys.exit(0)
#                 }
#                 volatility_gen * v_gen = itr->second.first;
                        v_gen = itr.first  # ###################################################
#                 float advance_pr = itr->second.second;
                        advance_pr = itr.second  # ################################################
#                 if (!v_gen->is_initialized()){
                        if not v_gen.is_initialized():
#                     v_gen->initialize(instance_count);
                            v_gen.initialize(instance_count)
#                 }
#                 /*
#                 if ( ((float) rand() / (float) RAND_MAX) < advance_pr ){
                        if float(random.uniform(0, 1)) < advance_pr:
#                     v_gen->advance();
                            v_gen.advance()
#                 }
#                 */
#                 int skip = ceil(1.0 / advance_pr);
                        skip = math.ceil(1.0 / advance_pr)
#                 if (((q_template._instantiationCount)%skip)==0){
                        if q_template._instantiationCount % skip == 0:
#                     cerr << "Instantiated queries=" << q_template._instantiationCount << "\n";
                            print(f'Instantiated queries={q_template._instantiationCount}')
#                     v_gen->advance();
                            v_gen.advance()
#                 }
#                 id = v_gen->next_rand_index();
                        id0 = v_gen.next_rand_index()
#             } else {
                    else:
#                 double r_value = model::generate_random(_distribution_type, instance_count);
                        r_value = Model.generate_random(self._distribution_type, instance_count)
#                 id = round(r_value * instance_count);
                        id0 = round(r_value * instance_count)
#             }

#             id = (id>=instance_count) ? (instance_count-1) : id;
                    if id0 >= instance_count:
                        id0 = instance_count - 1
#             result.append("<");
                    result += '<'
#             result.append(mdl._namespace_map.replace(_resource_type));
                    result += mdl._namespace_map.replace(self._resource_type)
#             result.append(boost::lexical_cast<string>(id));
                    result += str(id0)
#             result.append(">");
                    result += '>'
#             return result;
                    return result
#         } else {
                else:
#             string result = "";
                    result = ''
#             vector<string> * restricted_instances = mdl._type_map.get_instances(mdl._namespace_map.replace(_resource_type), mdl._namespace_map.replace(*_type_restriction));
                    restricted_instances = mdl._type_map.get_instances(mdl._namespace_map.replace(self._resource_type), mdl._namespace_map.replace(self._type_restriction))
#             instance_count = restricted_instances->size();
                    instance_count = len(restricted_instances)

#             unsigned int index = 0;
                    index = 0

#             if (_distribution_type==DISTRIBUTION_TYPES::DYNAMIC){
                    if self._distribution_type == DistributionTypes.DYNAMIC:
#                 map<string, pair<volatility_gen*, float> >::const_iterator itr = q_template._volatility_table.find(_dynamic_model_name);
                        try:
                            itr = q_template._volatility_table.find(self._dynamic_model_name)
#                 if (itr==q_template._volatility_table.end()){
                        except KeyError:
#                     cerr << "[mapping_m_t::generate()]\tDynamic model " << _dynamic_model_name << " does not exist..." << "\n";
                            print(f'[mapping_m_t::generate()]\tDynamic model {self._dynamic_model_name} does not exist...')
#                     exit(0);
                            sys.exit(0)
#                 }
#                 volatility_gen * v_gen = itr->second.first;
                        v_gen = itr.first  # ############################################
#                 float advance_pr = itr->second.second;
                        advance_pr = itr.second  # ###########################################
#                 if (!v_gen->is_initialized()){
                        if not v_gen.is_initialized():
#                     v_gen->initialize(instance_count);
                            v_gen.initialize(instance_count)
#                 }
#                 /*
#                 if ( ((float) rand() / (float) RAND_MAX) < advance_pr ){
#                     v_gen->advance();
#                 }
#                 */
#                 int skip = ceil(1.0 / advance_pr);
                        skip = math.ceil(1.0 / advance_pr)
#                 if (((q_template._instantiationCount)%skip)==0){
                        if q_template._instantiationCount % skip == 0:
#                     cerr << "Instantiated queries=" << q_template._instantiationCount << "\n";
                            print(f'Instantiated queries={q_template._instantiationCount}')
#                     v_gen->advance();
                            v_gen.advance()
#                 }

#                 index = v_gen->next_rand_index();
                        index = v_gen.next_rand_index()
#             } else {
                    else:
#                 double r_value = model::generate_random(_distribution_type, instance_count);
                        r_value = Model.generate_random(self._distribution_type, instance_count)
#                 index = round(r_value * instance_count);
                        index = round(r_value * instance_count)
#             }

#             index = (index>=instance_count) ? (instance_count-1) : index;
                    if index >= instance_count:
                        index = instance_count - 1
#             result.append("<");
                    result += '<'
#             result.append((*restricted_instances)[index]);
                    result += restricted_instances[index]
#             result.append(">");
                    result += '>'
#             delete restricted_instances;
#             return result;
                    return result
#         }
#     }
# }
        pass  # end of

    @staticmethod
# mapping_m_t * mapping_m_t::parse (const string & line){
    def parse(line):
        """

        Args:

        Returns:

        """
#     //cout<<"Parsing "<<line<<"\n";
#     string var_name;
        var_name = ''
#     bool is_literal_type = true;
        is_literal_type = True
#     LITERAL_TYPES::enum_t literal_type = LITERAL_TYPES::UNDEFINED;
        literal_type = LiteralTypes.UNDEFINED
#     string resource_type = "";
        resource_type = ''
#     bool type_restriction_exists = false;
        type_restriction_exists = False
#     string type_restriction = "";
        type_restriction = ''
#     DISTRIBUTION_TYPES::enum_t distribution_type = DISTRIBUTION_TYPES::UNDEFINED;
        distribution_type = DistributionTypes.UNDEFINED
#     string range_min = "";
        range_min = ''
#     string range_max = "";
        range_max = ''
#     string distribution_name = "";
        distribution_name = ''

#     stringstream parser(line);
        parser = parse_line(line)
#     int index = 0;
        index = 0
#     string token;
#     while (parser>>token){
        for token in parser:
#         switch (index){
#             case 0:{
            if index == 0:
#                 if (token.compare("#mapping")!=0){
                if token != '#mapping':
#                     cerr<<"[mapping_m_t::parse()]\tExpecting #mapping..."<<"\n";
                    print('[mapping_m_t::parse()]\tExpecting #mapping...')
#                     exit(0);
                    sys.exit(0)
#                 }
#                 break;
#             }
#             case 1:{
            elif index == 1:
    #                 var_name = token;
                var_name = token
#                 break;
#             }
#             case 2:{
            elif index == 2:
    #                 if (token.compare("integer")==0 || token.compare("INTEGER")==0){
                if token.find('integer') >= 0 or token.find('INTEGER') >= 0:
#                     literal_type = LITERAL_TYPES::INTEGER;
                    literal_type = LiteralTypes.INTEGER
#                 } else if (token.compare("string")==0 || token.compare("STRING")==0){
                elif token.find('string') >= 0 or token.find('STRING') >= 0:
#                     literal_type = LITERAL_TYPES::STRING;
                    literal_type = LiteralTypes.STRING
#                 } else if (token.compare("name")==0 || token.compare("NAME")==0){
                elif token.find('name') >= 0 or token.find('NAME') >= 0:
#                     literal_type = LITERAL_TYPES::NAME;
                    literal_type = LiteralTypes.NAME
#                 } else if (token.compare("date")==0 || token.compare("date")==0){
                elif token.find('date') >= 0 or token.find('DATE') >= 0:
#                     literal_type = LITERAL_TYPES::STRING;
                    literal_type = LiteralTypes.STRING
#                 } else {
                else:
#                     is_literal_type = false;
                    is_literal_type = False
#                     int pos = token.find_first_of('@');
                    pos = token.find('@')
#                     if (pos==string::npos){
                    if pos < 0:
#                         resource_type = token;
                        resource_type = token
#                     } else {
                    else:
#                         resource_type = token.substr(0, pos);
                        resource_type = token[:pos]
#                         type_restriction = token.substr(pos+1);
                        type_restriction = token[pos+1:]
#                         type_restriction_exists = true;
                        type_restriction_exists = True
#                     }
#                 }
#                 break;
#             }
#             case 3:{
            elif index == 3:
    #                 if (token.compare("uniform")==0 || token.compare("UNIFORM")==0){
                if token.find('uniform') >= 0 or token.find('UNIFORM') >= 0:
#                     distribution_type = DISTRIBUTION_TYPES::UNIFORM;
                    distribution_type = DistributionTypes.UNIFORM
#                 } else if (token.compare("normal")==0 || token.compare("NORMAL")==0){
                elif token.find('normal') >= 0 or token.find('NORMAL') >= 0:
#                     distribution_type = DISTRIBUTION_TYPES::NORMAL;
                    distribution_type = DistributionTypes.NORMAL
#                 } else if (token.compare("zipfian")==0 || token.compare("ZIPFIAN")==0){
                elif token.find('zipfian') >= 0 or token.find('ZIPFIAN') >= 0:
#                     distribution_type = DISTRIBUTION_TYPES::ZIPFIAN;
                    distribution_type = DistributionTypes.ZIPFIAN
#                 } else if (boost::starts_with(token, "dynamic#")==0 || boost::starts_with(token, "DYNAMIC#")==0){
                elif token.startswith('dynamic#') or token.startswith('DYNAMIC#'):
#                     distribution_type = DISTRIBUTION_TYPES::DYNAMIC;
                    distribution_type = DistributionTypes.DYNAMIC
#                     distribution_name = token.substr(8);
                    distribution_name = token[:8]  # #####################################
#                 }
#                 break;
#             }
#             case 4:{
            elif index == 4:
    #                 range_min = token;
                range_min = token
#                 break;
#             }
#             case 5:{
            elif index == 5:
        #                 range_max = token;
                range_max = token
#                 break;
#             }
#         }
#         index++;
            index += 1
#     }
#     if (is_literal_type){
        if is_literal_type:
#         if (index==3){
            if index == 3:
#             return new mapping_m_t(var_name, literal_type);
                return MappingMT(var_name=var_name, literal_type=literal_type)
#         } else if (index==4){
            elif index == 4:
#             return new mapping_m_t(var_name, literal_type, distribution_type);
                return MappingMT(var_name=var_name, literal_type=literal_type, distribution_type=distribution_type)
#         } else if (index==6){
            elif index == 6:
#             return new mapping_m_t(var_name, literal_type, distribution_type, range_min, range_max);
                return MappingMT(var_name=var_name, literal_type=literal_type, distribution_type=distribution_type, range_min=range_min, range_max=range_max)
#         }
#     } else {
        else:
#         if (type_restriction_exists){
            if type_restriction_exists:
#             if (index==3){
                if index == 3:
#                 return new mapping_m_t(var_name, resource_type, type_restriction);
                    return MappingMT(var_name=var_name, resource_type=resource_type, type_restriction=type_restriction)
#             } else if (index==4){
                elif index == 4:
#                 if (distribution_type==DISTRIBUTION_TYPES::DYNAMIC){
                    if distribution_type == DistributionTypes.DYNAMIC:
#                     mapping_m_t * result = new mapping_m_t(var_name, resource_type, type_restriction, distribution_type);
                        result = MappingMT(var_name=var_name, resource_type=resource_type, type_restriction=type_restriction, distribution_type=distribution_type)
#                     result->_dynamic_model_name = distribution_name;
                        result._dynamic_model_name = distribution_name
#                     return result;
                        return result
#                 } else {
                    else:
#                     return new mapping_m_t(var_name, resource_type, type_restriction, distribution_type);
                        return MappingMT(var_name=var_name, resource_type=resource_type, type_restriction=type_restriction, distribution_type=distribution_type)
#                 }
#             }
#         } else {
            else:
#             if (index==3){
                if index == 3:
#                 return new mapping_m_t(var_name, resource_type);
                    return MappingMT(var_name=var_name, resource_type=resource_type)
#             } else if (index==4){
                elif index == 4:
#                 if (distribution_type==DISTRIBUTION_TYPES::DYNAMIC){
                    if distribution_type == DistributionTypes.DYNAMIC:
#                     mapping_m_t * result = new mapping_m_t(var_name, resource_type, distribution_type);
                        result = MappingMT(var_name=var_name, resource_type=resource_type, distribution_type=distribution_type)
#                     result->_dynamic_model_name = distribution_name;
                        result._dynamic_model_name = distribution_name
#                     return result;
                        return result
#                 } else {
                    else:
#                     return new mapping_m_t(var_name, resource_type, distribution_type);
                        return MappingMT(var_name=var_name, resource_type=resource_type, distribution_type=distribution_type)
#                 }
#             }
#         }
#     }
#     cerr<<"[mapping_m_t::parse()]\tIncompatible arguments..."<<"\n";
        print('[mapping_m_t::parse()]\tIncompatible arguments...')
#     exit(0);
        sys.exit(0)
# }
        pass  # end of parse
    pass  # end of MappingMT


# struct operation_m_t {
class OperationMT:
    """

    Attributes:

    """
#     operation_m_t (const string & target_variable, const string & source_variable, OPERATION_TYPES::enum_t operation, int operand);
#     operation_m_t (const operation_m_t & rhs);

#     string compute (const map<string, string> & value_mappings); // @suppress("Invalid template argument")

#     static operation_m_t * parse (const string & line);
# };

# operation_m_t::operation_m_t (const string & target_variable, const string & source_variable, OPERATION_TYPES::enum_t operation, int operand){
    def __init__(self, target_variable=None, source_variable=None, operation=None, operand=None, rhs=None):
        """

        Args:

        Returns:

        """
#     string                      _target_variable;
        self._target_variable = ''
#     string                      _source_variable;
        self._source_variable = ''
#     OPERATION_TYPES::enum_t     _operation;
        self._operation = None
#     int                         _operand;
        self._operand = 0

        if rhs is None:
#     _target_variable = target_variable;
            self._target_variable = target_variable
#     _source_variable = source_variable;
            self._source_variable = source_variable
#     _operation = operation;
            self._operation = operation
#     _operand = operand;
            self._operand = operand
# }
#
# operation_m_t::operation_m_t (const operation_m_t & rhs){
        if rhs is not None:
#     _target_variable = rhs._target_variable;
            self._target_variable = rhs._target_variable
#     _source_variable = rhs._source_variable;
            self._source_variable = rhs._source_variable
#     _operation = rhs._operation;
            self._operation = rhs._operation
#     _operand = rhs._operand;
            self._operand = rhs._operand
# }
        pass  # end of __init__

# string operation_m_t::compute (const map<string, string> & value_mappings){
    def compute(self, value_mappings):
        """

        Args:

        Returns:

        """
#     int value = boost::lexical_cast<int>(value_mappings.find(_source_variable)->second);
        value = int(value_mappings[self._source_variable])
#     if (_operation==OPERATION_TYPES::ADDITION){
        if self._operation == OperationTypes.ADDITION:
#         value = value + _operand;
            value += self._operand
#     } else if (_operation==OPERATION_TYPES::MULTIPLICATION){
        elif self._operation == OperationTypes.MULTIPLICATION:
#         value = value * _operand;
            value *= self._operand
#     } else if (_operation==OPERATION_TYPES::MOD){
        elif self._operation == OperationTypes.MOD:
#         value = value % _operand;
            value %= self._operand
#     }
#     return boost::lexical_cast<string>(value);
        return str(value)
# }
        pass  # end of compute

    @staticmethod
# operation_m_t * operation_m_t::parse (const string & line){
    def parse(line):
        """

        Args:

        Returns:

        """
#     string target_variable = "", source_variable = "";
        target_variable = ''
        source_variable = ''
#     OPERATION_TYPES::enum_t operation = OPERATION_TYPES::UNDEFINED;
        operation = OperationTypes.UNDEFINED
#     int operand = -1;
        operand = -1
#     stringstream parser(line);
        parser = parse_line(line)
#     int index = 0;
        index = 0
#     string token;
#     while (parser>>token){
        for token in parser:
#         switch (index){
#             case 0:{
            if index == 0:
#                 if (token.compare("#operation")!=0){
                if token != '#operation':
#                     cerr<<"[operation_m_t::parse()]\tExpecting #operation..."<<"\n";
                    print('[operation_m_t::parse()]\tExpecting #operation...')
#                     exit(0);
                    sys.exit(0)
#                 }
#                 break;
#             }
#             case 1:{
            elif index == 1:
#                 target_variable = token;
                target_variable = token
#                 break;
#             }
#             case 2:{
            elif index == 2:
#                 source_variable = token;
                source_variable = token
#                 break;
#             }
#             case 3:{
            elif index == 3:
#                 if (token.compare("+")==0){
                if token.find('+') >= 0:
#                     operation = OPERATION_TYPES::ADDITION;
                    operation = OperationTypes.ADDITION
#                 } else if (token.compare("*")==0){
                elif token.find('*') >= 0:
#                     operation = OPERATION_TYPES::MULTIPLICATION;
                    operation = OperationTypes.MULTIPLICATION
#                 } else if (token.compare("mod")==0){
                elif token.find('mod') >= 0:
#                     operation = OPERATION_TYPES::MOD;
                    operation = OperationTypes.MOD
#                 }
#                 break;
#             }
#             case 4:{
            elif index == 4:
#                 operand = boost::lexical_cast<int>(token);
                operand = int(token)
#                 break;
#             }
#         }
#         index++;
            index += 1
#     }
#     if (index==5){
        if index == 5:
#         return new operation_m_t (target_variable, source_variable, operation, operand);
            return OperationMT(target_variable=target_variable, source_variable=source_variable, operation=operation, operand=operand)
#     } else {
        else:
#         cerr<<"[operation_m_t::parse()]\tIncompatible number of arguments..."<<"\n";
            print('[operation_m_t::parse()]\tIncompatible number of arguments...')
#         exit(0);
            sys.exit(0)
#     }
# }
        pass  # end of parse
    pass  # end of OperationMT


# struct query_template_m_t {
class QueryTemplateMT:
    """

    Attributes:

    """

#     query_template_m_t(const model * mdl);
#     query_template_m_t(const query_template_m_t & rhs);
#     ~query_template_m_t();

#     void instantiate (unsigned int query_count, unsigned int recurrence, vector<string> & result_array);

#     void parse (const string filename);
#     void parse_str (const string & content);
# };

# query_template_m_t::query_template_m_t(const model * mdl){
    def __init__(self, mdl=None, rhs=None):
        """

        Args:

        Returns:

        """
#     const model *                              _mdl;
        self._mdl = None
#     vector<mapping_m_t*>                        _variable_mapping_array;
        self._variable_mapping_array = []
#     vector<operation_m_t*>                      _operation_array;
        self._operation_array = []
#     vector<string>                              _template_lines;
        self._template_lines = []
    # map<string, pair<volatility_gen*, float> >  _volatility_table; // @suppress("Invalid template argument")
        self._volatility_table = {}
#     int                                         _instantiationCount;#     const model *                              _mdl;
        self._instantiationCount = 0

        if rhs is None:
#     _mdl = mdl;
            self._mdl = mdl
#     _instantiationCount = 0;
            self._instantiationCount = 0
# }

# query_template_m_t::query_template_m_t(const query_template_m_t & rhs){ // @suppress("Class members should be properly initialized")
        if rhs is not None:
#     for (vector<mapping_m_t*>::const_iterator itr=rhs._variable_mapping_array.cbegin(); itr!=rhs._variable_mapping_array.cend(); itr++){
            for itr in rhs._variable_mapping_array:
#         mapping_m_t * mapping = *itr;
                mapping = itr
#         _variable_mapping_array.push_back(new mapping_m_t(*mapping));
                self._variable_mapping_array.append(mapping)
#     }
#     for (vector<operation_m_t*>::const_iterator itr=rhs._operation_array.cbegin(); itr!=rhs._operation_array.cend(); itr++){
            for itr in rhs._operation_array:
#         operation_m_t * operation = *itr;
                operation = itr
#         _operation_array.push_back(new operation_m_t(*operation));
                self._operation_array.append(OperationMT(rhs=operation))
#     }
#     for (map<string, pair<volatility_gen*, float> >::const_iterator itr=rhs._volatility_table.begin(); itr!=rhs._volatility_table.end(); itr++){ // @suppress("Type cannot be resolved") // @suppress("Method cannot be resolved")
            for itr_key, itr_value in rhs._volatility_table.items():
#         string name = itr->first; // @suppress("Field cannot be resolved")
                name = itr_key
#         volatility_gen * v_gen = itr->second.first;
                v_gen = itr_value[0]
#         float advance_pr = itr->second.second;
                advance_pr = itr_value[1]
#         _volatility_table.insert(pair<string, pair<volatility_gen*,float> >(name, pair<volatility_gen*,float>(new volatility_gen(*v_gen), advance_pr)));
                self._volatility_table[name] = (VolatilityGen(v_gen), advance_pr)
#     }
#     _template_lines.insert(_template_lines.end(), rhs._template_lines.cbegin(), rhs._template_lines.cend());
            self._template_lines += rhs._template_lines
#     _instantiationCount = rhs._instantiationCount;
            self._instantiationCount = rhs._instantiationCount
# }
        pass  # end of query_template_m_t::query_template_m_t

# query_template_m_t::~query_template_m_t(){
#     for (vector<mapping_m_t*>::iterator itr=_variable_mapping_array.begin(); itr!=_variable_mapping_array.end(); itr++){
#         mapping_m_t * mapping = *itr;
#         delete mapping;
#     }
#     for (vector<operation_m_t*>::iterator itr=_operation_array.begin(); itr!=_operation_array.end(); itr++){
#         operation_m_t * operation = *itr;
#         delete operation;
#     }
#     for (map<string, pair<volatility_gen*, float> >::iterator itr=_volatility_table.begin(); itr!=_volatility_table.end(); itr++){
#         // volatility_gen * v_gen = itr->second.first;
#         //cout << "Deleting " << itr->first << " with address " << v_gen << "\n";
#         // delete v_gen;
#     }
# }

# void query_template_m_t::instantiate (unsigned int query_count, unsigned int recurrence, vector<string> & result_array){
    def instantiate(self, query_count, recurrence):
        """

        Args:

        Returns:

        """
        result_array = []
#     map<string, vector<string> > sample_map;
        sample_map = {}
#     unsigned int sample_count = (int)((float)query_count/(float)recurrence)+1;
        sample_count = int(float(query_count)/float(recurrence)) + 1
#     for (unsigned i=0; i<sample_count; i++){
        for i in range(sample_count):
#         for (vector<mapping_m_t*>::const_iterator itr=_variable_mapping_array.cbegin(); itr!=_variable_mapping_array.cend(); itr++){
            for itr in self._variable_mapping_array:
#             mapping_m_t * mapping = *itr;
                mapping = itr
#             if (mapping->_distribution_type!=DISTRIBUTION_TYPES::DYNAMIC){
                if mapping._distribution_type != DistributionTypes.DYNAMIC:
#                 if (sample_map.find(mapping->_var_name)==sample_map.end()){
                    try:
                        xxx = sample_map[mapping._var_name]
                    except KeyError:
#                     sample_map.insert(pair<string, vector<string> >(mapping->_var_name, vector<string>()));
                        sample_map[mapping._var_name] = []
#                 }
#                 sample_map[mapping->_var_name].push_back(mapping->generate(*_mdl, *this));
                    sample_map[mapping._var_name].append(mapping.generate(self._mdl, self))
#             }
#         }
#     }

#     for (unsigned i=0; i<query_count; i++){
        for i in range(query_count):
#         string query = "";
            query = ''
#         map<string, string> value_map;
            value_map = {}
#         ///////////////////////////////////////////////////////////////////////////////////////////
#         // Instead of populating value_map using mapping_m_t
#         //  just sample values from sample_map...
#         ///////////////////////////////////////////////////////////////////////////////////////////

#         for (vector<mapping_m_t*>::const_iterator itr=_variable_mapping_array.cbegin(); itr!=_variable_mapping_array.cend(); itr++){
            for itr in self._variable_mapping_array:
#             mapping_m_t * mapping = *itr;
                mapping = itr
#             if (mapping->_distribution_type==DISTRIBUTION_TYPES::DYNAMIC){
                if mapping._distribution_type == DistributionTypes.DYNAMIC:
#                 value_map.insert(pair<string, string>(mapping->_var_name, mapping->generate(*_mdl, *this)));
                    value_map[mapping._var_name] = mapping.generate(self._mdl, self)
#             } else {
                else:
#                 vector<string> samples = sample_map[mapping->_var_name];
                    samples = sample_map[mapping._var_name]
#                 value_map.insert(pair<string, string>(mapping->_var_name, samples[rand()%samples.size()]));
                    index = random.randint(0, len(samples)-1)
                    value_map[mapping._var_name] = samples[index]
#             }
#         }

#         /*
#         for (map<string, vector<string> >::iterator itr=sample_map.begin(); itr!=sample_map.end(); itr++){
#             vector<string> samples = itr->second;
#             value_map.insert(pair<string, string>(itr->first, samples[rand()%samples.size()]));
#         }
#         */

#         for (vector<operation_m_t*>::const_iterator itr=_operation_array.cbegin(); itr!=_operation_array.cend(); itr++){
            for itr in self._operation_array:
#             operation_m_t * operation = *itr;
                operation = itr
#             string value = operation->compute(value_map);
                value = operation.compute(value_map)
#             value_map.insert(pair<string, string>(operation->_target_variable, value));
                value_map[operation._target_variable] = value
#         }
#         for (vector<string>::const_iterator itr=_template_lines.cbegin(); itr!=_template_lines.cend(); itr++){
            for itr in self._template_lines:
#             string line = *itr;
                line = itr
#             float probability = 1.01;
                probability = 1.01

#             string line_cpy = line;
                line_cpy = line
#             boost::trim(line_cpy);
                line_cpy.strip()
#             if (line_cpy[0]=='['){
                if line_cpy[0] == '[':
#                 int begin_pos = line.find_first_of("[");
                    begin_pos = line.find('[')
#                 int end_pos = line.find_first_of("]");
                    end_pos = line.find(']')
#                 string line_prefix = line.substr(0, begin_pos);
                    line_prefix = line[:begin_pos]
#                 string probability_str = line.substr(begin_pos+1, end_pos-begin_pos-1);
                    probability_str = line[begin_pos+1: end_pos]
#                 string line_suffix = line.substr(end_pos+1);
                    line_suffix = line.substr[end_pos+1:]
#                 probability = boost::lexical_cast<float>(probability_str);
                    probability = float(probability_str)
#                 line = line_prefix;
                    line = line_prefix
#                 line.append(line_suffix);
                    line += line_suffix
#             }

#             /// Now you randomly generate a number, and check if it satisfies the probability.
#             /// If not you do not include this triple pattern in the query...
#             float random_f = ((float) rand()) / ((float) RAND_MAX);
                random_f = (float(random.uniform(0, 1)))  # / ((float) RAND_MAX)
#             if (random_f > probability){
                if random_f > probability:
#                 continue;
                    continue
#             }

#             string modified_line = "";
                modified_line = ''

#             // FIXME :: You need to implement a much better version of replace_all()...
#             string token;
#             int token_count = 0;
                token_count = 0
#             stringstream tokenizer (line);
                xxx = line.replace('\t', ' ')
                xxx = xxx.replace('  ', ' ')
                xxx = xxx.replace('  ', ' ')
                tokenizer = xxx.split(' ')
#             while (tokenizer>>token){
                for token in tokenizer:
#                 token_count++;
                    token_count += 1
#             }

#             int counter = 0;
                counter = 0
#             stringstream tokenizer2 (line);
                xxx = line.replace('\n', ' ')
                xxx = xxx.replace('\t', ' ')
                xxx = xxx.replace('  ', ' ')
                xxx = xxx.replace('  ', ' ')
                xxx.strip()
                tokenizer2 = xxx.split(' ')
#             while (tokenizer2>>token){
                for token in tokenizer2:
#                 if (token.find(":")!=string::npos){
                    if token.find(':') >= 0:
#                     modified_line.append("<");
                        modified_line += '<'
#                     modified_line.append(_mdl->_namespace_map.replace(token));
                        modified_line += self._mdl._namespace_map.replace(token)
#                     modified_line.append(">");
                        modified_line += '>'
#                 } else {
                    else:
#                     modified_line.append(token);
                        modified_line += token
#                 }
#                 counter++;
                    counter += 1
#                 if (counter<token_count){
                    if counter < token_count:
#                     modified_line.append(" ");
                        modified_line += ' '
#                 }
#             }

#             while (modified_line.find('%', 0)!=string::npos){
                while modified_line.find('%%') >= 0:
#                 int begin = modified_line.find('%', 0);
                    begin = modified_line.find('%%', 0)
#                 int end = modified_line.find('%', begin+1);
                    end = modified_line.find('%%', begin+1)
#                 string var_name = modified_line.substr(begin+1, end-begin-1);
                    var_name = modified_line[begin+1: end]
#                 modified_line = modified_line.replace(begin, end-begin+1, value_map[var_name]);
                    modified_line = modified_line.replace(f'%%{var_name}%%', value_map[var_name])
#             }
#             if (modified_line[modified_line.size()-1]=='.'){
                if modified_line[-1] == '.':
#                 query.append("\t");
                    query += '\t'
#             }
#             query.append(modified_line);
                query += modified_line
#             query.append("\n");
                query += '\n'
#         }
#         result_array.push_back(query);
            result_array.append(query)
#         _instantiationCount++;
            self._instantiationCount += 1
#     }
        return result_array
# }
        pass  # end of instantiate

# void query_template_m_t::parse (const string filename){
    def parse(self, filename):
        """

        Args:

        Returns:

        """
#     ifstream fis (filename);
        with (open(filename, 'r') as fis):
            lines = fis.readlines()
#     string line;
#     while (fis.good() && !fis.eof()){
#         getline(fis, line);
            for line in lines:
#         if (fis.good() && !fis.eof()){
#             if (boost::starts_with(line, "#mapping")){
                if line.startswith('#mapping'):
#                 mapping_m_t * mapping = mapping_m_t::parse(line);
                    mapping = MappingMT.parse(line)
#                 _variable_mapping_array.push_back(mapping);
                    self._variable_mapping_array.append(mapping)
#             } else if (boost::starts_with(line, "#operation")){
                elif line.startswith('#operation'):
#                 operation_m_t * operation = operation_m_t::parse(line);
                    operation = OperationMT.parse(line)
#                 _operation_array.push_back(operation);
                    self._operation_array.append(operation)
#             } else if (boost::starts_with(line, "#dynamic")){
                elif line.startswith('#dynamic'):
#                 string name = "";
                    name = ''
#                 float advance_pr = -1.0;
                    advance_pr = -1.0
#                 volatility_gen * v_gen = volatility_gen::parse(line, name, advance_pr);
                    v_gen = VolatilityGen.parse(line, name, advance_pr)
#                 if (_volatility_table.find(name)==_volatility_table.end()){
                    try:
                        xxx = self._volatility_table[name]
                        print('[query_template_m_t::parse]\tIgnoring duplicate #dynamic entry...')
                    except KeyError:
#                     _volatility_table.insert(pair<string, pair<volatility_gen*, float> >(name, pair<volatility_gen*, float>(v_gen, advance_pr)));
                        self._volatility_table[name] = (v_gen, advance_pr)
#                 } else {
#                     cerr << "[query_template_m_t::parse]\tIgnoring duplicate #dynamic entry..." << "\n";
#                 }
#             } else {
                else:
#                 _template_lines.push_back(line);
                    self._template_lines.append(line)
#             }
#         }
#     }
#     fis.close();
# }
        pass  # end of parse

# void query_template_m_t::parse_str (const string & content){
    def parse_str(self, content):
        """

        Args:

        Returns:

        """
#     stringstream sstream (content);
        lines = content.split('\n')
#     string line;
#     while (getline(sstream, line, '\n')){
        for line in lines:
#         if (boost::starts_with(line, "#mapping")){
            if line.startswith('#mapping'):
#             mapping_m_t * mapping = mapping_m_t::parse(line);
                mapping = MappingMT.parse(line)
#             _variable_mapping_array.push_back(mapping);
                self._variable_mapping_array.append(mapping)
#         } else if (boost::starts_with(line, "#operation")){
            elif line.startswith('#operation'):
#             operation_m_t * operation = operation_m_t::parse(line);
                operation = OperationMT.parse(line)
#             _operation_array.push_back(operation);
                self._operation_array.append(operation)
#         } else if (boost::starts_with(line, "#dynamic")){
            elif line.startswith('#dynamic'):
#                 string name = "";
                name = ''
#                 float advance_pr = -1.0;
                advance_pr = -1.0
#                 volatility_gen * v_gen = volatility_gen::parse(line, name, advance_pr);
                v_gen = None
                try:
                    v_gen = VolatilityGen.parse(line, name, advance_pr)
                    print('[query_template_m_t::parse]\tIgnoring duplicate #dynamic entry...')
#                 if (_volatility_table.find(name)==_volatility_table.end()){
                except KeyError:
#                     _volatility_table.insert(pair<string, pair<volatility_gen*, float> >(name, pair<volatility_gen*, float>(v_gen, advance_pr)));
                    self._volatility_table[name] = (v_gen, advance_pr)
#                 } else {
#                     cerr << "[query_template_m_t::parse]\tIgnoring duplicate #dynamic entry..." << "\n";
#                 }
#         } else {
            else:
#             _template_lines.push_back(line);
                self._template_lines.append(line)
#         }
#     }
# }
        pass  # end of parse_str
    pass  # end of QueryTemplateMT


# struct statistics_m_t {
class StatisticsMT:
    """

    Attributes:

    """
#     void init (const model * mdl, const string & predicate, const string & subject_type, const string & object_type);
#     statistics_m_t (const model * mdl, const string & predicate, const string & subject_type, const string & object_type);
#     statistics_m_t (const model * mdl, const string & predicate, const string & subject_type, const string & object_type, const string * subject_type_restriction, const string * object_type_restriction);
#     statistics_m_t (const statistics_m_t & rhs);
#     ~statistics_m_t();

#     void collect (const model * mdl, const string & subject, const string & predicate, const string & object);
#     void report () const;

#     static statistics_m_t * parse (const model * mdl, const string & line);
# };

# void statistics_m_t::init (const model * mdl, const string & predicate, const string & subject_type, const string & object_type){
    def init(self, mdl, predicate, subject_type, object_type):
        """

        Args:

        Returns:

        """
#     _predicate = predicate;
#     _subject_type = subject_type;
#     _object_type = object_type;
#     _subject_type_restriction = NULL;
#     _object_type_restriction = NULL;
        self._predicate = predicate
        self._subject_type = subject_type
        self._object_type = object_type
        self._subject_type_restriction = None
        self._object_type_restriction = None

#     _left_count = mdl->_id_cursor_map.find(subject_type)->second;
        self._left_count = mdl._id_cursor_map[subject_type]
#     _left_statistics = new int [_left_count];
        self._left_statistics = []
#     for (unsigned int i=0; i<_left_count; i++){
        for i in range(self._left_count):
#         _left_statistics[i] = 0;
            self._left_statistics.append(0)
#     }
#     _right_count = mdl->_id_cursor_map.find(object_type)->second;
        self._right_count = mdl._id_cursor_map[object_type]
#     _right_statistics = new int [_right_count];
        self._right_statistics = []
#     for (unsigned int i=0; i<_right_count; i++){
        for i in range(self._right_count):
#         _right_statistics[i] = 0;
            self._right_statistics.append(0)
#     }
# }
        pass  # end of init

# statistics_m_t::statistics_m_t (const model * mdl, const string & predicate, const string & subject_type, const string & object_type){
    def __init__(self, mdl=None, predicate=None, subject_type=None, object_type=None, subject_type_restriction=None, object_type_restriction=None, rhs=None):
        """

        Args:

        Returns:

        """
#     string              _predicate;
#     string              _subject_type;
#     string              _object_type;
#     string *            _subject_type_restriction;
#     string *            _object_type_restriction;
#     unsigned int        _left_count;
#     int *               _left_statistics;
#     unsigned int        _right_count;
#     int *               _right_statistics;
        self._predicate = ''
        self._subject_type = ''
        self._object_type = ''
        self._subject_type_restriction = ''
        self._object_type_restriction = ''
        self._left_count = 0
        self._left_statistics = 0
        self._right_count = 0
        self._right_statistics = 0

        if mdl is not None and predicate is not None and subject_type is not None and object_type is not None:
#     init (mdl, predicate, subject_type, object_type);
            self.init(mdl, predicate, subject_type, object_type)
# }

# statistics_m_t::statistics_m_t (const model * mdl, const string & predicate, const string & subject_type, const string & object_type, const string * subject_type_restriction, const string * object_type_restriction){
#     init (mdl, predicate, subject_type, object_type);
        if subject_type_restriction is not None and object_type_restriction is not None:
#     if (subject_type_restriction!=NULL){
#         _subject_type_restriction = new string (*subject_type_restriction);
            self._subject_type_restriction = subject_type_restriction
#     }
#     if (object_type_restriction!=NULL){
#         _object_type_restriction = new string (*object_type_restriction);
            self._object_type_restriction = object_type_restriction
#     }
# }

# statistics_m_t::statistics_m_t (const statistics_m_t & rhs){
        if rhs is not None:
#     _predicate = rhs._predicate;
            self._predicate = rhs._predicate
#     _subject_type = rhs._subject_type;
            self._subject_type = rhs._subject_type
#     _object_type = rhs._object_type;
            self._object_type = rhs._object_type
#     if (rhs._subject_type_restriction!=NULL){
            if rhs._subject_type_restriction is not None:
#         _subject_type_restriction = new string (*rhs._subject_type_restriction);
                self._subject_type_restriction = rhs._subject_type_restriction
#     } else {
            else:
#         _subject_type_restriction = NULL;
                self._subject_type_restriction = None
#     }
#     if (rhs._object_type_restriction!=NULL){
            if rhs._object_type_restriction is not None:
#         _object_type_restriction = new string (*rhs._object_type_restriction);
                self._object_type_restriction = rhs._object_type_restriction
#     } else {
            else:
#         _object_type_restriction = NULL;
                self._object_type_restriction = None
#     }
#     _left_count = rhs._left_count;
            self._left_count = rhs._left_count
#     _left_statistics = new int [_left_count];
            self._left_statistics = []
#     for (unsigned int i=0; i<_left_count; i++){
            for i in range(self._left_count):
#         _left_statistics[i] = rhs._left_statistics[i];
                self._left_statistics.append(rhs._left_statistics[i])
#     }
#     _right_count = rhs._right_count;
            self._right_count = rhs._right_count
#     _right_statistics = new int [_right_count];
            self._right_statistics = []
#     for (unsigned int i=0; i<_right_count; i++){
            for i in range(self._right_count):
#         _right_statistics[i] = rhs._right_statistics[i];
                self._right_statistics.append(rhs._right_statistics[i])
#     }
# }
        pass  # end of __init__

# statistics_m_t::~statistics_m_t(){
#     delete _subject_type_restriction;
#     delete _object_type_restriction;
#     delete [] _left_statistics;
#     delete [] _right_statistics;
# }

# void statistics_m_t::collect (const model * mdl, const string & subject, const string & predicate, const string & object){
    def collect(self, mdl, subject, predicate, object0):
        """

        Args:

        Returns:

        """
#     string p_uri = "";
        p_uri = ''
#     p_uri.append("<");
        p_uri += '<'
#     p_uri.append(mdl->_namespace_map.replace(_predicate));
        p_uri += mdl._namespace_map.replace(self._predicate)
#     p_uri.append(">");
        p_uri += '>'
#     if (p_uri.compare(predicate)==0){
        if p_uri == predicate:
#         string s_prefix = "", o_prefix = "";
            s_prefix = ''
            o_prefix = ''
#         s_prefix.append("<");
            s_prefix += '<'
#         s_prefix.append(mdl->_namespace_map.replace(_subject_type));
            s_prefix += mdl._namespace_map.replace(self._subject_type)
#         o_prefix.append("<");
            o_prefix += '<'
#         o_prefix.append(mdl->_namespace_map.replace(_object_type));
            o_prefix += mdl._namespace_map.replace(self._object_type)
#         if (boost::starts_with(subject, s_prefix) && boost::starts_with(object, o_prefix)){
            if subject.startswith(s_prefix) and object0.startswith(o_prefix):
#             // FIXME::Also check for type restrictions if they exist...
#             string s_instance = subject.substr(1, subject.size()-2);
                s_instance = subject[1:-2]
#             string o_instance = object.substr(1, object.size()-2);
                o_instance = object0[1:-2]
#             unsigned int s_id = boost::lexical_cast<unsigned int>(subject.substr(s_prefix.size(), (subject.size() - s_prefix.size() - 1)));
                s_id = int(subject[len(s_prefix), (len(subject) - len(s_prefix) - 1)])
#             unsigned int o_id = boost::lexical_cast<unsigned int>(object.substr(o_prefix.size(), (object.size() - o_prefix.size() - 1)));
                o_id = int(object0[len(o_prefix), (len(object0) - len(o_prefix) - 1)])
#             if ((_subject_type_restriction==NULL || mdl->_type_map.instanceof(s_instance, mdl->_namespace_map.replace(*_subject_type_restriction))) &&
#                 (_object_type_restriction==NULL || mdl->_type_map.instanceof(o_instance, mdl->_namespace_map.replace(*_object_type_restriction)))){
                if (self._subject_type_restriction is None or mdl._type_map.instanceof(s_instance, mdl._namespace_map.replace(self._subject_type_restriction))) and \
                   (self._object_type_restriction is None or mdl._type_map.instanceof(o_instance, mdl._namespace_map.replace(self._object_type_restriction))):
#                 _left_statistics[s_id] = _left_statistics[s_id] + 1;
                    self._left_statistics[s_id] += 1
#                 _right_statistics[o_id] = _right_statistics[o_id] + 1;
                    self._right_statistics[s_id] += 1
#             } else {
                else:
#                 if (!(_subject_type_restriction==NULL || mdl->_type_map.instanceof(s_instance, mdl->_namespace_map.replace(*_subject_type_restriction)))){
                    if not (self._subject_type_restriction is None or mdl._type_map.instanceof(s_instance, mdl._namespace_map.replace(self._subject_type_restriction))):
#                     _left_statistics[s_id] = -1;
                        self._left_statistics[s_id] = -1
#                 }
#                 if (!(_object_type_restriction==NULL || mdl->_type_map.instanceof(o_instance, mdl->_namespace_map.replace(*_object_type_restriction)))){
                    if not (self._object_type_restriction is None or mdl._type_map.instanceof(o_instance, mdl._namespace_map.replace(self._object_type_restriction))):
#                     _right_statistics[o_id] = -1;
                        self._right_statistics[o_id] = -1
#                 }
#             }
#         }
#     }
# }
        pass  # end of collect

# void statistics_m_t::report () const{
    def report(self):
        """

        Args:

        Returns:

        """
#     float left_cover = 0.0;
        left_cover = 0.0
#     unsigned int left_actual_count = 0;
        left_actual_count = 0
#     unsigned int left_min = numeric_limits<unsigned int>::max();
        left_min = 1000000
#     unsigned int left_max = numeric_limits<unsigned int>::min();
        left_max = -1000000
#     float left_mean = 0.0;
        left_mean = 0.0
#     vector<unsigned int> left_distribution;
        left_distribution = []

#     for (unsigned int i=0; i<_left_count; i++){
        for i in range(self._left_count):
#         if (_left_statistics[i]>0){
            if self._left_statistics[i] > 0:
#             left_cover += 1.0;
                left_cover += 1.0
#             if (_left_statistics[i]<left_min){
                if self._left_statistics[i] < left_min:
#                 left_min = _left_statistics[i];
                    left_min = self._left_statistics[i]
#             }
#             if (_left_statistics[i]>left_max){
                if self._left_statistics[i] > left_max:
#                 left_max = _left_statistics[i];
                    left_max = self._left_statistics[i]
#             }
#             left_mean += _left_statistics[i];
                left_mean += self._left_statistics[i]
#             left_distribution.push_back(_left_statistics[i]);
                left_distribution.append(self._left_statistics[i])
#         }
#         if (_left_statistics[i]>=0){
            if self._left_statistics[i] >= 0:
#             ++left_actual_count;
                left_actual_count += 1
#         }
#     }
#     left_mean = left_mean / left_cover;
        left_mean = left_mean / left_cover
#     left_cover = left_cover / ((float) left_actual_count);
        left_cover = left_cover / float(left_actual_count)
#     sort(left_distribution.begin(), left_distribution.end());
        left_distribution.sort()

#     float right_cover = 0.0;
        right_cover = 0.0
#     unsigned int right_actual_count = 0;
        right_actual_count = 0
#     unsigned int right_min = numeric_limits<unsigned int>::max();
        right_min = 1000000
#     unsigned int right_max = numeric_limits<unsigned int>::min();
        right_max = -1000000
#     float right_mean = 0.0;
        right_mean = 0.0
#     vector<unsigned int> right_distribution;
        right_distribution = []

#     for (unsigned int i=0; i<_right_count; i++){
        for i in range(self._right_count):
#         if (_right_statistics[i]>0){
            if self._left_statistics[i] > 0:
#             right_cover += 1.0;
                right_cover += 1.0
#             if (_right_statistics[i]<right_min){
                if self._right_statistics[i] < right_min:
#                 right_min = _right_statistics[i];
                    right_min = self._right_statistics[i]
#             }
#             if (_right_statistics[i]>right_max){
                if self._right_statistics[i] > right_max:
#                 right_max = _right_statistics[i];
                    right_max = self._right_statistics[i]
#             }
#             right_mean += _right_statistics[i];
                right_mean += self._right_statistics[i]
#             right_distribution.push_back(_right_statistics[i]);
                right_distribution.append(self._right_statistics[i])
#         }
#         if (_right_statistics[i]>=0){
            if self._right_statistics[i] >= 0:
#             ++right_actual_count;
                right_actual_count += 1
#         }
#     }
#     right_mean = right_mean / right_cover;
        right_mean = right_mean / right_cover
#     right_cover = right_cover / ((float) right_actual_count);
        right_cover = right_cover / float(right_actual_count)
#     sort(right_distribution.begin(), right_distribution.end());
        right_distribution.sort()

#     cout<<"Printing statistics..."<<"\n";
        print('Printing statistics...')
#     cout<<"\tPredicate:           "<<_predicate<<"\n";
        print(f'\tPredicate:           {self._predicate}')
#     cout<<"\tSubject-type:        "<<_subject_type<<"\n";
        print(f'\tSubject-type:           {self._subject_type}')
#     cout<<"\tObject-type:         "<<_object_type<<"\n";
        print(f'\tObject-type:           {self._object_type}')
#     if (_subject_type_restriction!=NULL){
        if self._subject_type_restriction is not None:
#         cout<<"\tSubject-restriction: "<<*_subject_type_restriction<<"\n";
            print(f'\tSubject-restriction: {self._subject_type_restriction}')
#     }
#     if (_object_type_restriction!=NULL){
        if self._object_type_restriction is not None:
#         cout<<"\tObject-restriction: "<<*_object_type_restriction<<"\n";
            print(f'\tObject-restriction: {self._subject_type_restriction}')
#     }

#     cout<<"\t\tSubject-statistics..."<<"\n";
        print('\t\tSubject-statistics...')
#     cout<<"\t\t\tCover:        "<<left_cover<<"\n";
        print(f'\t\t\tCover:        {left_cover}')
#     cout<<"\t\t\tRange:        "<<"["<<left_min<<"-"<<left_max<<"]"<<"\n";
        print(f'\t\t\tRange:        [{left_min}-{left_max}]')
#     cout<<"\t\t\tMean:         "<<left_mean<<"\n";
        print(f'\t\t\tMean:         {left_mean}')
#     cout<<"\t\t\tDistribution: ";
        print('\t\t\tDistribution: ')
#     for (vector<unsigned int>::reverse_iterator itr=left_distribution.rbegin(); itr!=left_distribution.rend(); itr++){
        for itr in left_distribution:
#         cout<<*itr<<" ";
            print(itr, ' ', end='')
#     }
#     cout<<"\n";
        print()

#     cout<<"\t\tObject-statistics..."<<"\n";
        print('\t\tObject-statistics...')
#     cout<<"\t\t\tCover:        "<<right_cover<<"\n";
        print(f'\t\t\tCover:        {right_cover}')
#     cout<<"\t\t\tRange:        "<<"["<<right_min<<"-"<<right_max<<"]"<<"\n";
        print(f'\t\t\tRange:        [{right_min}-{right_max}]')
#     cout<<"\t\t\tMean:         "<<right_mean<<"\n";
        print(f'\t\t\tMean:         {right_mean}')
#     cout<<"\t\t\tDistribution: ";
        print('\t\t\tDistribution: ')
#     for (vector<unsigned int>::reverse_iterator itr=right_distribution.rbegin(); itr!=right_distribution.rend(); itr++){
        for itr in right_distribution:
#         cout<<*itr<<" ";
            print(itr, ' ', end='')
#     }
#     cout<<"\n";
        print()
# }
        pass  # end of report
    @staticmethod
# statistics_m_t * statistics_m_t::parse (const model * mdl, const string & line){
    def parse(mdl, line):
        """

        Args:

        Returns:

        """
#     string subject_type ('');
        subject_type = ''
#     string predicate ("");
        predicate = ''
#     string object_type ("");
        object_type = ''
#     string * subject_type_restriction = NULL;
        subject_type_restriction = None
#     string * object_type_restriction = NULL;
        object_type_restriction = None

#     stringstream parser(line);
        parser = parse_line(line)
#     int index = 0;
        index = 0
#     string token;
#     while (parser>>token){
        for token in parser:
#         switch (index){
#             case 0: {
            if index == 0:
#                 if (token.compare("#statistics")!=0){
                if token.find('#statistics') < 0:
#                     cerr<<"[statistics_m_t::parse()]\tExpecting #statistics..."<<"\n";
                    print('[statistics_m_t::parse()]\tExpecting #statistics...')
#                     exit(0);
                    sys.exit(0)
#                 }
#                 break;
#             }
#             case 1: {
            elif index == 1:
#                 subject_type = token;
                subject_type = token
#                 break;
#             }
#             case 2: {
            elif index == 2:
#                 predicate = token;
                predicate = token
#                 break;
#             }
#             case 3: {
            elif index == 3:
#                 object_type = token;
                object_type = token
#                 break;
#             }
#             case 4: {
            elif index == 4:
#                 if (!boost::starts_with(token, "@")){
                if token.startswith('@'):
#                     cerr<<"[statistics_m_t::parse()]\tExpecting '@' qualifier before type restriction..."<<"\n";
                    print('[statistics_m_t::parse()]\tExpecting '@' qualifier before type restriction...')
#                     exit(0);
                    sys.exit(0)
#                 }
#                 if (token.compare("@null")!=0 && token.compare("@NULL")!=0){
                if token.find('@null') < 0 and token.find('@NULL') < 0:
#                     subject_type_restriction = new string(token.substr(1));
                    subject_type_restriction = str(token[1:])
#                 }
#                 break;
#             }
#             case 5: {
            elif index == 5:
#                 if (!boost::starts_with(token, "@")){
                if token.startswith('@'):
#                     cerr<<"[statistics_m_t::parse()]\tExpecting '@' qualifier before type restriction..."<<"\n";
                    print('[statistics_m_t::parse()]\tExpecting '@' qualifier before type restriction...')
#                     exit(0);
                    sys.exit(0)
#                 }
#                 if (token.compare("@null")!=0 && token.compare("@NULL")!=0){
                if token.find('@null') < 0 and token.find('@NULL') < 0:
#                     object_type_restriction = new string(token.substr(1));
                    object_type_restriction = str(token[1:])
#                 }
#                 break;
#             }
#         }
#         index++;
            index += 1
#     }
#     statistics_m_t * result = NULL;
        result = None
#     if (index==4){
        if index == 4:
#         result = new statistics_m_t (mdl, predicate, subject_type, object_type);
            result = StatisticsMT(mdl=mdl, predicate=predicate, subject_type=subject_type, object_type=object_type)
#     } else if (index==6){
        elif index == 6:
#         result = new statistics_m_t (mdl, predicate, subject_type, object_type, subject_type_restriction, object_type_restriction);
            result = StatisticsMT(mdl=mdl, predicate=predicate, subject_type=subject_type, object_type=object_type, subject_type_restriction=subject_type_restriction, object_type_restriction=object_type_restriction)
#     } else {
        else:
#         cerr<<"[statistics_m_t::parse()]\tExpecting 3 or 5 arguments..."<<"\n";
            print('[statistics_m_t::parse()]\tExpecting 3 or 5 arguments...')
#         exit(0);
            sys.exit(0)
#     }
#     delete subject_type_restriction;
#     delete object_type_restriction;
#     return result;
        return result
# }
        pass  # end of parse
    pass # end of StatisticsMT


class Model:
    """

    Attributes:

    """
    zipfian_cache = {}
# model::model(const char * filename){

    def __init__(self, filename: str):
        """

        Args:

        Returns:

        """
        # vector<resource_m_t*>       _resource_array;
        self._resource_array: list[ResourceMT] = []
        # vector<association_m_t*>    _association_array;
        self._association_array = []
        # vector<string>              _statistics_lines;
        self._statistics_lines = []
        # map<string, unsigned int>   _id_cursor_map; // @suppress("Invalid template argument")
        self._id_cursor_map: dict[str, int] = {}
        # namespace_map               _namespace_map;
        self._namespace_map: NamespaceMap = NamespaceMap()
        # type_map                    _type_map;
        self._type_map: TypeMap = TypeMap()
        #     srand (time(NULL));
        #     parse(filename);
        self.parse(filename)
    # }
        pass

# model::~model(){
#     for (vector<resource_m_t*>::iterator itr=_resource_array.begin(); itr!=_resource_array.end(); itr++){
#         delete *itr;
#     }
#     for (vector<association_m_t*>::iterator itr=_association_array.begin(); itr!=_association_array.end(); itr++){
#         delete *itr;
#     }
# }

# void model::generate (int scale_factor){
    def generate(self, scale_factor: int):
        """

        Args:

        Returns:

        """
#     boost::posix_time::ptime t1 (bpt::microsec_clock::universal_time());

#     for (int i=0; i<scale_factor; i++){
        for i in range(scale_factor):
#         for (vector<resource_m_t*>::iterator itr2=_resource_array.begin(); itr2!=_resource_array.end(); itr2++){
            for itr2 in self._resource_array:
#             resource_m_t * resource = *itr2;
                resource = itr2
#             if (i==0 || resource->_scalable){
                if i == 0 or resource._scalable:
#                 resource->generate(_namespace_map, _id_cursor_map);
                    resource.generate(self._namespace_map, self._id_cursor_map)
#             }
                    pass  # end of if
#         }
                pass  # end of for
#     }
            pass  # end of for

#     boost::posix_time::ptime t2 (bpt::microsec_clock::universal_time());

#     for (vector<association_m_t*>::iterator itr1=_association_array.begin(); itr1!=_association_array.end(); itr1++){
        for itr1 in self._association_array:
#         association_m_t * association = *itr1;
            association: AssociationMT = itr1
#         association->generate(_namespace_map, _type_map, _id_cursor_map);
            association.generate(self._namespace_map, self._type_map, self._id_cursor_map)
#     }
        pass  # end of for

#     boost::posix_time::ptime t3 (bpt::microsec_clock::universal_time());

#     for (vector<resource_m_t*>::iterator itr1=_resource_array.begin(); itr1!=_resource_array.end(); itr1++){
        for itr1 in self._resource_array:
#         resource_m_t * resource = *itr1;
            resource = itr1
#         resource->process_type_restrictions(_namespace_map, _type_map, _id_cursor_map);
            resource.process_type_restrictions(self._namespace_map, self._type_map, self._id_cursor_map)
#     }
        pass  # end of for

#     boost::posix_time::ptime t4 (bpt::microsec_clock::universal_time());

#     for (vector<association_m_t*>::iterator itr1=_association_array.begin(); itr1!=_association_array.end(); itr1++){
        for itr1 in self._association_array:
#         association_m_t * association = *itr1;
            association = itr1
#         association->process_type_restrictions(_namespace_map, _type_map, _id_cursor_map);
            association.process_type_restrictions(self._namespace_map, self._type_map, self._id_cursor_map)
#     }
        pass  # end of for

#     boost::posix_time::ptime t5 (bpt::microsec_clock::universal_time());

#     //cerr << "[t1--t2]" << " " << (t2-t1).total_microseconds() << "\n";
#     //cerr << "[t2--t3]" << " " << (t3-t2).total_microseconds() << "\n";
#     //cerr << "[t3--t4]" << " " << (t4-t3).total_microseconds() << "\n";
#     //cerr << "[t4--t5]" << " " << (t5-t4).total_microseconds() << "\n";
# }
        pass  # end of model::generate

# void model::compute_statistics (const vector<triple_st> & triples){
    def compute_statistics(self, triples):
        """

        Args:

        Returns:

        """
#     vector<statistics_m_t*> statistics_array;
        statistics_array = []
#     for (vector<string>::iterator itr=_statistics_lines.begin(); itr!=_statistics_lines.end(); itr++){
        for itr in self._statistics_lines:
#         statistics_m_t * statistics = statistics_m_t::parse(this, *itr);
            statistics = StatisticsMT.parse(self, itr)
#         statistics_array.push_back(statistics);
            statistics_array.append(statistics)
#     }
#     if (!statistics_array.empty()){
        if len(statistics_array) == 0:
#         for (vector<triple_st>::const_iterator itr1=triples.begin(); itr1!=triples.end(); itr1++){
            for itr1 in triples:
#             for (vector<statistics_m_t*>::iterator itr2=statistics_array.begin(); itr2!=statistics_array.end(); itr2++){
                for itr2 in statistics_array:
#                 statistics_m_t * statistics = *itr2;
                    statistics = itr2
#                 statistics->collect(this, itr1->_subject, itr1->_predicate, itr1->_object);
                    statistics.collect(self, itr1._subject, itr1._predicate, itr1._object)
#             }
#         }
#         for (vector<statistics_m_t*>::iterator itr=statistics_array.begin(); itr!=statistics_array.end(); itr++){
            for itr in statistics_array:
#             statistics_m_t * statistics = *itr;
                statistics = itr
#             statistics->report();
                statistics.report()
#             delete statistics;
#         }
#     }
# }
        pass

# void model::parse (const char * filename){
    def parse(self, filename):
        """

        Args:

        Returns:

        """
#     stack<pair<short,void*> > object_stack;
        object_stack = []
#     ifstream fis (filename);
        with open(filename, 'r') as fis:
            lines = fis.readlines()
#     string line;
#     while (fis.good() && !fis.eof()){
            for line in lines:
                line = line.replace('\n', '')
#         getline(fis, line);
#         if (boost::starts_with(line, "//")){
                if line.startswith('//'):
#             // Ignore the line, it is a comment...
                    pass
#         } else if (boost::starts_with(line, "#namespace")){
                elif line.startswith('#namespace'):
#             namespace_m_t * namespace_declaration = namespace_m_t::parse(line);
                    namespace_declaration = NamespaceMT.parse(line)
#             _namespace_map.insert(*namespace_declaration);
                    self._namespace_map.insert(namespace_declaration=namespace_declaration)
#             delete namespace_declaration;
#         } else if (boost::starts_with(line, "#predicate")){
                elif line.startswith('#predicate'):
#             predicate_m_t * predicate = predicate_m_t::parse(line);
                    predicate = PredicateMT.parse(line)
#             object_stack.push(pair<short,void*>(2, (void*)predicate));
                    object_stack.append((2, predicate))
#         } else if (boost::starts_with(line, "<pgroup>")){
                elif line.startswith('<pgroup>'):
#             predicate_group_m_t * predicate_group = predicate_group_m_t::parse(line);
                    predicate_group = PredicateGroupMT.parse(line)
#             object_stack.push(pair<short,void*>(1, (void*)predicate_group));
                    object_stack.append((1, predicate_group))
#         } else if (boost::starts_with(line, "</pgroup>")){
                elif line.startswith('</pgroup>'):
#             vector<predicate_m_t*> array;
                    array = []
#             while (object_stack.top().first!=1){
                    while object_stack[-1][0] != 1:
#                 array.push_back((predicate_m_t * ) object_stack.top().second);
                        array.append(object_stack[-1][1])
#                 object_stack.pop();
                        object_stack = object_stack[:-1]
#             }
#             predicate_group_m_t * predicate_group = (predicate_group_m_t *) object_stack.top().second;
                    predicate_group = object_stack[-1][1]
#             predicate_group->_predicate_array.insert(predicate_group->_predicate_array.begin(), array.begin(), array.end());
                    predicate_group._predicate_array = array + predicate_group._predicate_array
#         } else if (boost::starts_with(line, "<type")){
                elif line.startswith('<type'):
#             resource_m_t * resource = resource_m_t::parse(line);
                    resource = ResourceMT.parse(line)
#             object_stack.push(pair<short,void*>(0, (void*)resource));
                    object_stack.append((0, resource))
#         } else if (boost::starts_with(line, "</type>")){
                elif line.startswith('</type>'):
#             vector<predicate_group_m_t*> array;
                    array = []
#             while (object_stack.top().first!=0){
                    while len(object_stack) > 0 and object_stack[-1][0] != 0:
#                 array.push_back((predicate_group_m_t * ) object_stack.top().second);
                        array.append(object_stack[-1][1])
#                 object_stack.pop();
                        object_stack = object_stack[:-1]
#             }
#             resource_m_t * resource = (resource_m_t *) object_stack.top().second;
                    resource = object_stack[-1][1]
#             resource->_predicate_group_array.insert(resource->_predicate_group_array.begin(), array.begin(), array.end());
                    resource._predicate_group_array = array + resource._predicate_group_array
#         } else if (boost::starts_with(line, "#association")){
                elif line.startswith('#association'):
#             association_m_t * association = association_m_t::parse(_id_cursor_map, line);
                    association = AssociationMT.parse(line)
#             _association_array.push_back(association);
                    self._association_array.append(association)
#         } else if (boost::starts_with(line, "#statistics")){
                elif line.startswith('#statistics'):
#             _statistics_lines.push_back(line);
                    self._statistics_lines.append(line)
#         }
#     }
#     fis.close();
#     while (!object_stack.empty()){
        while len(object_stack) > 0:
#         _resource_array.push_back((resource_m_t *) object_stack.top().second);
            self._resource_array.append(object_stack[-1][1])
#         object_stack.pop();
            object_stack = object_stack[:-1]
#     }
# }
        pass

    @staticmethod
# string model::generate_literal (LITERAL_TYPES::enum_t literal_type, DISTRIBUTION_TYPES::enum_t distribution_type, const string & range_min, const string & range_max){
    def generate_literal(literal_type, distribution_type, range_min, range_max):
        """

        Args:

        Returns:

        """
#     string literal = "";
        literal = ''
#     switch (literal_type){
#         case LITERAL_TYPES::INTEGER:{
        if literal_type == LiteralTypes.INTEGER:
#             int min_value = boost::lexical_cast<int>(range_min);
            min_value = int(range_min)
#             int max_value = boost::lexical_cast<int>(range_max);
            max_value = int(range_max)
#             int interval = max_value - min_value;
            interval = max_value - min_value
#             double r_value = model::generate_random(distribution_type);
            r_value = Model.generate_random(distribution_type)
#             int offset = round(r_value * interval);
            offset = round(r_value * interval)
#             offset = (offset<0) ? 0 : offset;
            if offset < 0:
                offset = 0
#             offset = (offset>interval) ? interval : offset;
            if offset > interval:
                offset = interval
#             literal.append(boost::lexical_cast<string>(min_value + offset));
            literal += str(min_value + offset)
#             break;
#         }
#         case LITERAL_TYPES::STRING:{
        elif literal_type == LiteralTypes.STRING:
#             // pair<unsigned int, unsigned int> range = dictionary::get_instance()->get_interval(DICTIONARY_TYPES::ENGLISH_WORDS, range_min, range_max);
#             // int interval = range.second - range.first - 1;
#             // double r_value = model::generate_random(distribution_type);
#             // int offset = round(r_value * interval);
#             // offset = (offset<0) ? 0 : offset;
#             // offset = (offset>interval) ? interval : offset;
#             // literal.append(*(dictionary::get_instance()->get_word(DICTIONARY_TYPES::ENGLISH_WORDS, range.first + offset)));
#             // // Keep appending a few more words from the dictionary...
#             // unsigned int wc = rand() % MAX_LITERAL_WORDS;
#             // unsigned int dict_size = dictionary::get_instance()->word_count(DICTIONARY_TYPES::ENGLISH_WORDS);
#             // for (unsigned int index=0; index<wc; index++){
#             //     literal.append(" ");
#             //     literal.append(*(dictionary::get_instance()->get_word(DICTIONARY_TYPES::ENGLISH_WORDS, rand()%dict_size)));
#             // }
#             break;
            pass
#         }
#         case LITERAL_TYPES::NAME:{
        elif literal_type == LiteralTypes.NAME:
#             // pair<unsigned int, unsigned int> range = dictionary::get_instance()->get_interval(DICTIONARY_TYPES::FIRST_NAMES, range_min, range_max);
#             // int interval = range.second - range.first - 1;
#             // double r_value = model::generate_random(distribution_type);
#             // int offset = round(r_value * interval);
#             // offset = (offset<0) ? 0 : offset;
#             // offset = (offset>interval) ? interval : offset;
#             // literal.append(*(dictionary::get_instance()->get_word(DICTIONARY_TYPES::FIRST_NAMES, range.first + offset)));
#             break;
            pass
#         }
#         case LITERAL_TYPES::DATE:{
        elif literal_type == LiteralTypes.DATE:
#             boost::posix_time::ptime min_time, max_time, gen_time;
#             locale format (locale::classic(),new boost::posix_time::time_input_facet("%Y-%m-%d"));
#             istringstream min_iss (range_min);
#             min_iss.imbue(format);
#             min_iss>>min_time;
#             istringstream max_iss (range_max);
#             max_iss.imbue(format);
#             max_iss>>max_time;
#             boost::posix_time::time_duration range (max_time - min_time);
#             long interval = range.total_seconds();
#             double r_value = model::generate_random(distribution_type);
#             long offset = round(r_value * interval);
#             offset = (offset<0) ? 0 : offset;
#             offset = (offset>interval) ? interval : offset;
#             gen_time = min_time + boost::gregorian::days(offset/(24*3600));
#             boost::gregorian::date gen_date = gen_time.date();
#             literal.append(boost::gregorian::to_iso_extended_string(gen_date));
#             break;
            pass
#         }
#     }
#     return literal;
        return literal
# }
        pass

    @staticmethod
# double model::generate_random (DISTRIBUTION_TYPES::enum_t distribution_type, int item_count){
    def generate_random(distribution_type, item_count=None):
        """

        Args:

        Returns:

        """

        # double model::generate_zipfian (int item_count){
        def generate_zipfian(item_count: int) -> float:
            """

            Args:

            Returns:

            """
            #     vector<double> * intervals = NULL;
            intervals = []

            #     if (zipfian_cache.find(item_count)==zipfian_cache.end()){
            try:
                intervals = Model.zipfian_cache[item_count]
            except KeyError:
            #         intervals = new vector<double>();
                intervals = []
            #         double offset = 0.0;
                offset = 0.0
            #         for (int i=1; i<=item_count; i++){
                for i in range(1, item_count):
            #             offset += 1.0 / ((double) (i));
                    offset += 1.0 / float(i)
            #             intervals->push_back(offset);
                    intervals.append(offset)
            #         }
            #         double scale_factor = 1.0 / offset;
                scale_factor = 1.0 / offset
            #         for (int cursor=0; cursor<item_count; cursor++){
                for cursor in range(item_count):
            #             (*intervals)[cursor] = (*intervals)[cursor] * scale_factor;
                    intervals[cursor] *= scale_factor
            #         }
            #         zipfian_cache.insert(pair<int, vector<double>*>(item_count, intervals));
                    zipfian_cache[item_count] = intervals
            #     } else {

            #         intervals = zipfian_cache[item_count];
            #     }
            #
            #     double random_value = ((double) rand()) / ((double) RAND_MAX);
            random_value = random.uniform(0, 1)
            #     vector<double>::iterator pivot = lower_bound(intervals->begin(), intervals->end(), random_value);
            pivot = bisect.bisect_left(intervals, random_value)
            #     double result = (pivot - intervals->begin()) * (1.0 / ((double) item_count));
            result = (pivot - intervals[0]) * (1.0 / float(item_count))
            #     return result;
            return result
            # }
            pass

#     double result = 0.0;
        result = 0.0
#     switch (distribution_type){
#         case DISTRIBUTION_TYPES::UNIFORM:{
        if distribution_type == DistributionTypes.UNIFORM:
#             result = ((double) rand()) / ((double) RAND_MAX);
            result = random.uniform(0, 1)
#             break;
#         }
#         case DISTRIBUTION_TYPES::NORMAL:{
        elif distribution_type == DistributionTypes.NORMAL:
#             result = BOOST_NORMAL_DIST_GEN();
            xxx = np.random.normal(0.5, 0.5/3.0, 1)
            result = xxx[0]
#             break;
#         }
#         case DISTRIBUTION_TYPES::ZIPFIAN:{
        elif distribution_type == DistributionTypes.ZIPFIAN:
#             result = generate_zipfian(item_count);
            result = generate_zipfian(item_count)
#             break;
#         }
#         case DISTRIBUTION_TYPES::UNDEFINED:
        elif distribution_type == DistributionTypes.UNDEFINED:
#         default:{
#             break;
            pass
#         }
#     }
#     result = (result<0.0) ? 0.0:result;
        if result < 0.0:
            result = 0.0
#     result = (result>1.0) ? 1.0:result;
        if result > 1.0:
            result = 1.0
#     return result;
        return result
# }
        pass

    @staticmethod
# void model::clear_zipfian_cache (){
    def clear_zipfian_cache():
        """

        Args:

        Returns:

        """
#     for (map<int, vector<double>*>::iterator itr=zipfian_cache.begin(); itr!=zipfian_cache.end(); itr++){
        for itr_key, itr_value in Model.zipfian_cache.items():
#         vector<double> * value = itr->second;
            value = itr_value
#         delete value;
#         itr->second = NULL;
#     }
#     zipfian_cache.clear();
        Model.zipfian_cache = {}
# }
        pass

# void model::load (const char * filename){
    def load(self, filename):
        """

        Args:

        Returns:

        """
#     _id_cursor_map.clear();
        self._id_cursor_map = {}
#     _type_map.clear();
        self._type_map = {}
#
#     ifstream fis (filename);
        with open(filename, 'r') as fis:
#     int lCount = -1;
            l_count = -1
#     string line, token;
            lines = fis.readlines()
#     getline(fis, line);
            line = lines[0]
            lines = lines[1:]
#     lCount = boost::lexical_cast<int>(line);
            l_count = int(line)
#     for (int i=0; i<lCount; i++){
            for i in range(l_count):
#         getline(fis, line);
                line = lines[0]
                lines = lines[1:]
#         stringstream parser(line);
                parser = line.split(' ')
#         parser>>token;
                token = parser[0]
#         string key = token;
                key = token
#         parser>>token;
                token = parser[1]
#         unsigned int value = boost::lexical_cast<unsigned int>(token);
                value = int(token)
#         _id_cursor_map.insert(pair<string, unsigned int>(key, value));
                self._id_cursor_map[key] = value
#     }

#     getline(fis, line);
            line = lines[0]
            lines = lines[1:]
#     lCount = boost::lexical_cast<int>(line);
            l_count = int(line)
#     for (int i=0; i<lCount; i++){
            for i in range(l_count):
#         getline(fis, line);
                line = lines[0]
                lines = lines[1:]
#         stringstream parser(line);
                parser = line.split(' ')
#         parser>>token;
                token = parser[0]
#         string type = token;
                type = token
#         while (parser>>token){
                for token in parser[1:]:
#             _type_map.insert(token, type);
                    self._type_map[token] = type
#         }
#     }
#
#     // You do not need to load namespaces...
#     // They come automatically from the model file...
#     /*
#     getline(fis, line);
#     lCount = boost::lexical_cast<int>(line);
#     for (int i=0; i<lCount; i++){
#         string key, value;
#         getline(fis, line);
#         stringstream parser(line);
#         parser>>key;
#         parser>>value;
#         _namespace_map.insert(key, value);
#     }
#     */
#
#     fis.close();
# }
        pass  # end of load

# void model::save (const char * filename) const{
    def save(self, filename):
        """

        Args:

        Returns:

        """
#     ofstream fos (filename);
        lines = []
        # with open(filename, 'w') as fos:

#     fos<<_id_cursor_map.size()<<"\n";
#             fos.write(str(len(self._id_cursor_map)))
        lines.append(str(len(self._id_cursor_map))+'\n')
#     for (map<string, unsigned int>::const_iterator itr1=_id_cursor_map.cbegin(); itr1!=_id_cursor_map.cend(); itr1++){
        for itr1_key, itr1_value in self._id_cursor_map.items():
#         fos<<itr1->first<<" "<<boost::lexical_cast<string>(itr1->second)<<"\n";
#             fos.write(itr1.first, itr1.second)
            lines.append(str(itr1_key)+' '+str(itr1_value)+'\n')
#     }
            pass  # end of for
        with open(filename, 'w') as fos:
            fos.writelines(lines)
    #     vector<string> lines;
#     _type_map.to_str(lines);
            lines = self._type_map.to_str()
#     fos<<lines.size()<<"\n";
            fos.write(str(len(lines))+'\n')
#     for (vector<string>::iterator itr1=lines.begin(); itr1!=lines.end(); itr1++){
            for itr1 in lines:
#         fos<<*itr1<<"\n";
                fos.write(itr1+'\n')
#     }
        pass  # end of save

#     // You do not need to save namespaces...
#     // They come automatically from the model file...
#     /*
#     lines.clear();
#     _namespace_map.to_str(lines);
#     fos<<lines.size()<<"\n";
#     for (vector<string>::iterator itr1=lines.begin(); itr1!=lines.end(); itr1++){
#         fos<<*itr1<<"\n";
#     }
#     */

#     fos.close();
# }
    pass


# int main(int argc, const char* argv[]) {
def main():
    """
    Main entry point to start the whole program.

    Args:
        arguments are given from command line.
    Returns:

    """
    argv = sys.argv  # arguments are given from command line.
#     dictionary * dict = dictionary::get_instance();
    dict0: Dictionary = Dictionary.get_instance()  # create an instance of Dictionary class
#     string xxx = argv[1];  ////
    argc = len(argv)  # number of arguments. argv[0] holds the program name.
#     if ( (argc==2 || argc==4 || argc==5 || argc>=6) && argv[1][0]=='-'){
    if (argc == 2 or argc == 4 or argc == 5 or argc >= 6) and argv[1][0] == '-':
#         dict->init("/usr/share/dict/words", "./data/files/firstnames.txt", "./data/files/lastnames.txt");
        dict0.init('/usr/share/dict/words', './data/files/firstnames.txt', './data/files/lastnames.txt')  # initialize the dictionary instance
#         const char * model_filename = argv[2];
        model_filename = argv[2]  # file name of the model
#         model cur_model (model_filename);
        cur_model = Model(model_filename)  # read the model
#         // statistics stat (cur_model);
#         if (argc==4 && argv[1][0]=='-' && argv[1][1]=='d'){
        if argc == 4 and argv[1] == '-d':  # '-d' for preparing saved.txt
#             unsigned int scale_factor = boost::lexical_cast<unsigned int>(string(argv[3]));
            scale_factor = int(argv[3])
#             cur_model.generate(scale_factor);
            cur_model.generate(scale_factor)
#             cur_model.save("saved.txt");
            cur_model.save('./data/saved.txt')
#             //statistics stat (&cur_model, triples);
#             // dictionary::destroy_instance();
#             return 0;
            return 0
#         } else if (argc==5 && argv[1][0]=='-' && argv[1][1]=='q'){
        elif argc == 5 and argv[1] == '-q':  # '-q' for generating the instances of queries from query templates
#             cur_model.load("saved.txt");
            cur_model.load('./data/saved.txt')
#             unsigned int query_count = boost::lexical_cast<unsigned int>(string(argv[(argc-2)]));
            query_count = int(str(argv[(argc-2)]))
#             unsigned int recurrence_factor = boost::lexical_cast<unsigned int>(string(argv[(argc-1)]));
            recurrence_factor = int(str(argv[(argc-1)]))
#             vector<string> workload;
            workload: list[str] = []
#             string line, qTemplateStr = "";
            qTemplateStr = ''
#             while (getline(cin, line)){
            while True:
                line = input()
#                 if (boost::starts_with(line, "#end")){
                if line.startswith('#end'):
#                     query_template_m_t q_template (&cur_model);
                    q_template = QueryTemplateMT(cur_model)
#                     q_template.parse_str(qTemplateStr);
                    q_template.parse_str(qTemplateStr)
#                     q_template.instantiate(query_count, recurrence_factor, workload);
                    workload = q_template.instantiate(query_count, recurrence_factor)
#                     qTemplateStr = "";
                    qTemplateStr = ''
                    break
#                 } else {
                else:
#                     qTemplateStr.append(line);
                    qTemplateStr += line
#                     qTemplateStr.append("\n");
                    qTemplateStr += '\n'
#                 }
#             }
#
#             // obtain a time-based seed:
#             //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
#             //shuffle (workload.begin(), workload.end(), std::default_random_engine(seed));
#
#             for (int qid=0; qid<workload.size(); qid++){
            for query in workload:
#                 cout<<workload[qid];
                print(query)
#             }
#             dictionary::destroy_instance();
#             return 0;
            return 0
#         } else if (argc>=6 && argv[1][0]=='-' && argv[1][1]=='q'){
        elif argc >= 6 and argv[1] == '-q':  # '-q' for generating the instances of queries from query templates
#             cur_model.load("saved.txt");
            cur_model.load('./data/saved.txt')  # read model schema
#             unsigned int query_count = boost::lexical_cast<unsigned int>(string(argv[(argc-2)]));
            query_count = int(str(argv[(argc-2)]))  # number of queries to be generated
#             unsigned int recurrence_factor = boost::lexical_cast<unsigned int>(string(argv[(argc-1)]));
            recurrence_factor = int(str(argv[(argc-1)]))
#             vector<string> workload;
            workload = []
#             for (int template_id=3; template_id<(argc-2); template_id++){
            for template_id in range(3, argc-2):  # file names for query templates
#                 const char * query_filename = argv[template_id];
                query_filename = argv[template_id]
#                 query_template_m_t q_template (&cur_model);
                q_template = QueryTemplateMT(mdl=cur_model)  # create an instance of QueryTemplateMT class
#                 q_template.parse(query_filename);
                q_template.parse(query_filename)  # read the query template file
#                 q_template.instantiate(query_count, recurrence_factor, workload);
                workload = q_template.instantiate(query_count, recurrence_factor)  # receive the generated queries
#             }
#
#             /// obtain a time-based seed:
#             //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
#             //shuffle (workload.begin(), workload.end(), std::default_random_engine(seed));
#
#             for (int qid=0; qid<workload.size(); qid++){
                for query in workload:
#                 cout<<workload[qid];
                    print(query)
#             }
#             // dictionary::destroy_instance();
#             return 0;
            return 0
#         } else if (argc==6 && argv[1][0]=='-' && argv[1][1]=='s'){
        elif argc == 6 and argv[1] == '-s':  # '-s' generates query templates
#             cur_model.load("saved.txt");
            cur_model.load('./data/saved.txt')
#             vector<triple_st> triple_array = triple_st::parse_file(argv[3]);
            triple_array = TripleST.parse_file(argv[3])  # create an instance of TripleST class
#             int maxQSize = boost::lexical_cast<int>(argv[4]);
            maxQSize = int(argv[4])  # max number of triples in a query
#             int qCount = boost::lexical_cast<int>(argv[5]);
            qCount = int(argv[5])  # number of queries to be generated
#             // statistics stat (&cur_model, triple_array, maxQSize, qCount, 1, false, false);
            Statistics(cur_model, triple_array, maxQSize, qCount, 1, False, False)  # generate and print query templates
#             // dictionary::destroy_instance();
#             return 0;
            return 0
#         } else if (argc==7 && argv[1][0]=='-' && argv[1][1]=='s'){
        elif argc == 7 and argv[1] == '-s':  # '-s' with an additional argument of constCount
#             cur_model.load("saved.txt");
            cur_model.load('./data/saved.txt')
#             vector<triple_st> triple_array = triple_st::parse_file(argv[3]);
            triple_array = TripleST.parse_file(argv[3])
#             int maxQSize = boost::lexical_cast<int>(argv[4]);
            maxQSize = int(argv[4])
#             int qCount = boost::lexical_cast<int>(argv[5]);
            qCount = int(argv[5])
#             int constCount = boost::lexical_cast<int>(argv[6]);
            constCount = int(argv[6])
#             // statistics stat (&cur_model, triple_array, maxQSize, qCount, constCount, false, false);
            Statistics(cur_model, triple_array, maxQSize, qCount, constCount, False, False)
#             // dictionary::destroy_instance();
#             return 0;
            return 0
#         } else if (argc==8 && argv[1][0]=='-' && argv[1][1]=='s'){
        elif argc == 8 and argv[1] == '-s':  # '-s' with an additional argument of constant-join-vertex-allowed?
#             cur_model.load("saved.txt");
            cur_model.load('./data/saved.txt')
#             vector<triple_st> triple_array = triple_st::parse_file(argv[3]);
            triple_array: TripleST = TripleST.parse_file(argv[3])
#             int maxQSize = boost::lexical_cast<int>(argv[4]);
            maxQSize: int = int(argv[4])
#             int qCount = boost::lexical_cast<int>(argv[5]);
            qCount: int = int(argv[5])
#             int constCount = boost::lexical_cast<int>(argv[6]);
            constCount: int = int(argv[6])
#             // statistics stat (&cur_model, triple_array, maxQSize, qCount, constCount, argv[7][0]=='t', false);
            Statistics(cur_model, triple_array, maxQSize, qCount, constCount, argv[7][0]=='t', False)
#             // dictionary::destroy_instance();
#             return 0;
            return 0
#         } else if (argc==9 && argv[1][0]=='-' && argv[1][1]=='s'){
        elif argc == 9 and argv[1] == '-s':  # '-s' with an additional argument of duplicate-edges-allowed? give 't'.
#             cur_model.load("saved.txt");
            cur_model.load('./data/saved.txt')
#             vector<triple_st> triple_array = triple_st::parse_file(argv[3]);
            triple_array = TripleST.parse_file(argv[3])
#             int maxQSize = boost::lexical_cast<int>(argv[4]);
            maxQSize = int(argv[4])
#             int qCount = boost::lexical_cast<int>(argv[5]);
            qCount = int(argv[5])
#             int constCount = boost::lexical_cast<int>(argv[6]);
            constCount = int(argv[6])
#             // statistics stat (&cur_model, triple_array, maxQSize, qCount, constCount, argv[7][0]=='t', argv[8][0]=='t');
            Statistics(cur_model, triple_array, maxQSize, qCount, constCount, argv[7][0] == 't', argv[8][0] == 't')
#             // dictionary::destroy_instance();
#             return 0;
            return 0
#         } else if (argc==2 && argv[1][0]=='-' && argv[1][1]=='x'){
        elif argc == 2 and argv[1] == '-x':  # for test.
#             // volatility_gen::test();
            VolatilityGen.test()
#             // dictionary::destroy_instance();
#             return 0;
            return 0
#         }
#         return 0;
        return 0
#     }
        pass  # end of if statement
    else:
#     std::cout<<"Usage:::\t./watdiv -d <model-file> <scale-factor>"<<"\n";
        print('Usage:::\t./watdiv -d <model-file> <scale-factor>')
#     std::cout<<"Usage:::\t./watdiv -q <model-file> <query-count> <recurrence-factor>"<<"\n";
        print('Usage:::\t./watdiv -q <model-file> <query-count> <recurrence-factor>')
#     std::cout<<"        \t./watdiv -q <model-file> <query-file> <query-count> <recurrence-factor>"<<"\n";
        print('        \t./watdiv -q <model-file> <query-file> <query-count> <recurrence-factor>')
#     std::cout<<"Usage:::\t./watdiv -s <model-file> <dataset-file> <max-query-size> <query-count>"<<"\n";
        print('Usage:::\t./watdiv -s <model-file> <dataset-file> <max-query-size> <query-count>')
#     std::cout<<"        \t./watdiv -s <model-file> <dataset-file> <max-query-size> <query-count> <constant-per-query-count>"<<"\n";
        print('        \t./watdiv -s <model-file> <dataset-file> <max-query-size> <query-count> <constant-per-query-count>')
#     std::cout<<"        \t./watdiv -s <model-file> <dataset-file> <max-query-size> <query-count> <constant-per-query-count> <constant-join-vertex-allowed?>"<<"\n";
        print('        \t./watdiv -s <model-file> <dataset-file> <max-query-size> <query-count> <constant-per-query-count> <constant-join-vertex-allowed?>')
#     //cout<<"Usage:::\t./watdiv -x"<<"\n";
#     //cout<<"        \t./watdiv -s <model-file> <dataset-file> <max-query-size> <query-count> <constant-per-query-count> <constant-join-vertex-allowed?> <duplicate-edges-allowed?>"<<"\n";
#     // dictionary::destroy_instance();
#     return 0;
        return 0
# }
    pass


if __name__ == '__main__':
    main()  # call the main entry point of the program
