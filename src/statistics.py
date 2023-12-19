"""
statistics.py
Python version of WatDiv/statistics
T. Masuda, 2023/12/1
The codes in this file have been moved to model.py
"""

import random
from enum import Enum, auto
# #include "statistics.h"
from src.model import DistributionTypes, LiteralTypes

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
#     string _vertex1;
#     string _vertex2;
#     string _edge;
#     DISTRIBUTION_TYPES::enum_t _right_distribution;
#     int _group_id;
#     double _pr_left;
#     double _pr_right;


# statistics_st::statistics_st (const string & vertex1, const string & vertex2, const string & edge, DISTRIBUTION_TYPES::enum_t right_distribution, int group_id){
    def __init__(self, vertex1, vertex2, edge, right_distribution, group_id):
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
        pass

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
# statistics::statistics(const model * mdl, const vector<triple_st> & triple_array, int maxQSize, int qCount, int constCount, bool constJoinVertexAllowed, bool dupEdgesAllowed){
    def __init__(self, mdl, triple_array, maxQSize, qCount, constCount, constJoinVertexAllowed, dupEdgesAllowed):
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

#     for (int i=0; i<qCount; ){
        for i in range(qCount):
#         int qSize = (rand() % maxQSize) + 1;
            qSize = int(random.uniform(0, maxQSize-1)) + 1
#         if (traverse_graph(qSize, constCount, constJoinVertexAllowed, dupEdgesAllowed, query_map)){
            if self.traverse_graph(qSize, constCount, constJoinVertexAllowed, dupEdgesAllowed, query_map):
#             i++;
                i += 1
#         }
#     }

#     /*
#     // Randomly select a subset of queries from each category...
#     cout << "Randomly selected query templates ..." << "\n";
        print("Randomly selected query templates ...")
#     map<pair<QUERY_STRUCTURE::enum_t, QUERY_CATEGORY::enum_t>, map<string, string> >::iterator itr1 = query_map.begin();
#     for (; itr1!=query_map.end(); itr1++){
        for itr1_key, itr1_value in self.query_map.items():
#         if (itr1->first.first != QUERY_STRUCTURE::UNDEFINED && itr1->first.second != QUERY_CATEGORY::UNDEFINED){
            if itr1_key[0] != QueryStructure.UNDEFINED and itr1_key[1] != QueryCategory.UNDEFINED:
#             cout << "CATEGORY ::: " << itr1->first.first << " - " << itr1->first.second << "\n";
                print(f'CATEGORY ::: {itr1_key[0]} - {itr1_key[1]}')
#             map<string, string> query_description_map = itr1->second;
                query_description_map = itr1_value
#             vector<string> keys;
                keys = []
#             for (map<string, string>::iterator itr2=query_description_map.begin(); itr2!=query_description_map.end(); itr2++){
                for itr2 in query_description_map:
#                 keys.push_back(itr2->first);
                    keys.append(itr2[0])
#             }
#             for (int k=0; k<5; k++){
                for k in range(5):
                    xxx = random.uniform(0, len(keys) - 1)
#                 string random_query = keys[rand()%keys.size()];
                    random_query = keys[xxx]
#                 string description = query_description_map[random_query];
                    description = query_description_map[random_query]
#                 cout << random_query << description << "\n";
                    print(f'{random_query} {description}')
#             }
#         }
#     }
            pass  # end of for
#     */
# }
        pass  # end of statistics::statistics

# statistics::~statistics(){
# }

# void statistics::extract_schema(const model & mdl){
    def extract_schema(self, mdl):
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
        pass  # end of

# void statistics::populate_graph(const vector<statistics_st> & tuples){
    def populate_graph(self, tuples):
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
#     for (map<string, set<statistics_st> >::const_iterator itr1=graph.begin(); itr1!=graph.end(); itr1++){
        for itr1_key, iter1_value in self.graph.items():
#         string vertex = itr1->first;
            vertex = itr1_key
#         int pos = string::npos;
            pos = vertex.find('@')
#         if ( (pos = vertex.find("@"))!=string::npos){
            if pos >= 0:
#             string base_vertex = vertex.substr(0, pos);
                base_vertex = vertex[0, pos]
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
#     vector<string> v_array;
        v_array = []
#     for (map<string, set<statistics_st> >::const_iterator itr1=graph.begin(); itr1!=graph.end(); itr1++){
        for itr1_key, itr1_value in self.graph.items():
#         v_array.push_back(itr1->first);
            v_array.append(itr1_key)
#     }

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
            potential_edges =[]
#         map<string, set<statistics_st> >::const_iterator fit = graph.find(v);
            fit = self.graph[v]
#         set<statistics_st> incident_edges = fit->second;
            incident_edges = fit
#         for (set<statistics_st>::const_iterator itr1=incident_edges.begin(); itr1!=incident_edges.end(); itr1++){
            for itr1 in incident_edges:
                found = False
                try:
                    xxx = traversed_edges[itr1]
                    found = True
                except KeyError:
                    pass
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
#         //cout << "\t" << "(" << next_edge._vertex1 << ", " << next_edge._edge << ", " << next_edge._vertex2 << ")" << "\n";
#         traversed_edges.insert(next_edge);

#         vector<string> potential_next_vertices;
#         if (graph.find(next_edge._vertex1)!=graph.end()){
#             potential_next_vertices.push_back(next_edge._vertex1);
#         }
#         if (graph.find(next_edge._vertex2)!=graph.end()){
#             potential_next_vertices.push_back(next_edge._vertex2);
#         }
#         v = potential_next_vertices[rand() % potential_next_vertices.size()];

#         itr_counter++;
#     } while (traversed_edges.size()<max_size && itr_counter<(max_size*20));

#     // Compute vertices in query graph...
#     set<string> variable_set;
#     string query_str = "";
#     int var_count = 0;
#     map<string, int> q_vertex_map;
#     map<string, set<string> > variable_map;
#     for (set<statistics_st>::iterator itr1=traversed_edges.begin(); itr1!=traversed_edges.end(); itr1++){
#         string var1 = "", var2 = "";
#         string v1_base = itr1->_vertex1, v2_base = itr1->_vertex2;
#         int pos = string::npos;
#         if ((pos=v1_base.find("@"))!=string::npos){
#             v1_base = v1_base.substr(0, pos);
#         }
#         if ((pos=v2_base.find("@"))!=string::npos){
#             v2_base = v2_base.substr(0, pos);
#         }
#         if (q_vertex_map.find(v1_base)==q_vertex_map.end()){
#             q_vertex_map.insert(pair<string, int>(v1_base, var_count));
#             var_count++;
#         }
#         var1.append("?v");
#         var1.append(boost::lexical_cast<string>(q_vertex_map[v1_base]));
#         query_str.append("\t");
#         query_str.append(var1);
#         query_str.append("\t");
#         query_str.append(itr1->_edge);
#         query_str.append("\t");
#         if (v2_base.compare("date") !=0 && v2_base.compare("integer") !=0 && v2_base.compare("name") !=0 && v2_base.compare("string")){
#             if (q_vertex_map.find(v2_base)==q_vertex_map.end()){
#                 q_vertex_map.insert(pair<string, int>(v2_base, var_count));
#                 var_count++;
#             }
#             var2.append("?v");
#             var2.append(boost::lexical_cast<string>(q_vertex_map[v2_base]));
#             query_str.append(var2);
#         } else {
#             var2.append("?v");
#             var2.append(boost::lexical_cast<string>(var_count));
#             query_str.append(var2);
#             var_count++;
#         }
#         query_str.append(" ");
#         query_str.append(".");
#         query_str.append("\n");
#         if (variable_map.find(var1)==variable_map.end()){
#             variable_map.insert(pair<string, set<string> >(var1, set<string>()));
#         }
#         variable_map[var1].insert(itr1->_vertex1);
#         if (variable_map.find(var2)==variable_map.end()){
#             variable_map.insert(pair<string, set<string> >(var2, set<string>()));
#         }
#         variable_map[var2].insert(itr1->_vertex2);
#         variable_set.insert(var1);
#         variable_set.insert(var2);
#     }

#     /////////////////////////////////////////////////////////////////////////////
#     // Put all traversed edges in a graph...
#     // When constructing the graph, be careful to ignore literals...
#     // For each join vertex (i.e., vertex with more than one incident edge),
#     //  compute join selectivity.
#     /////////////////////////////////////////////////////////////////////////////

#     var_count = 0;
#     map<string, set<statistics_st> > query_graph;
#     for (set<statistics_st>::iterator itr1 = traversed_edges.begin(); itr1!= traversed_edges.end(); itr1++){
#         string v1_base = itr1->_vertex1, v2_base = itr1->_vertex2;
#         int pos = string::npos;
#         if ((pos=v1_base.find("@"))!=string::npos){
#             v1_base = v1_base.substr(0, pos);
#         }
#         if ((pos=v2_base.find("@"))!=string::npos){
#             v2_base = v2_base.substr(0, pos);
#         }
#         if (query_graph.find(v1_base)==query_graph.end()){
#             query_graph.insert(pair<string, set<statistics_st> >(v1_base, set<statistics_st>()));
#         }
#         query_graph[v1_base].insert(*itr1);

#         if (v2_base.compare("date") !=0 && v2_base.compare("integer") !=0 && v2_base.compare("name") !=0 && v2_base.compare("string")){
#             if (query_graph.find(v2_base)==query_graph.end()){
#                 query_graph.insert(pair<string, set<statistics_st> >(v2_base, set<statistics_st>()));
#             }
#             query_graph[v2_base].insert(*itr1);
#         } else {
#             string tmp_id = "?v";
#             tmp_id.append(boost::lexical_cast<string>(var_count));
#             query_graph[tmp_id].insert(*itr1);
#             var_count++;
#         }
#     }

#     /*
#     QUERY_STRUCTURE::enum_t query_structure = QUERY_STRUCTURE::UNDEFINED;
#     QUERY_CATEGORY::enum_t query_category = QUERY_CATEGORY::UNDEFINED;
#     double min_join_selectivity = 1.0;
#     double max_join_selectivity = 1.0;
#     */

#     /*
#     vector<int> degree_array;
#     for (map<string, set<statistics_st> >::iterator itr1=query_graph.begin(); itr1!=query_graph.end(); itr1++){
#         const set<statistics_st> edge_list = itr1->second;
#         degree_array.push_back(edge_list.size());
#         if (edge_list.size()>=2){
#             set<int> correlation_set;
#             for (set<statistics_st>::const_iterator itr2=edge_list.begin(); itr2!=edge_list.end(); itr2++){
#                 if (correlation_set.find(itr2->_group_id)==correlation_set.end()){
#                     correlation_set.insert(itr2->_group_id);
#                     string entity = "";
#                     bool direction = true;
#                     if (itr2->_vertex1.find(itr1->first)!=string::npos){
#                         entity = itr2->_vertex1;
#                         direction = true;
#                     } else if (itr2->_vertex2.find(itr1->first)!=string::npos){
#                         entity = itr2->_vertex2;
#                         direction = false;
#                     } else {
#                         cout<<"Direction is wrong..."<<"\n";
#                         cout<<"Vertex "<<itr1->first<<"\n";
#                         cout<<"Edge-1 "<<itr2->_vertex1<<"\n";
#                         cout<<"Edge-2 "<<itr2->_vertex2<<"\n";
#                         exit(0);
#                     }
#                     double max_stats = 0.0;
#                     double min_stats = 1.0;
#                     for (int dist=0; dist<DISTRIBUTION_TYPES::UNDEFINED; dist++){
#                         string key = get_key(entity, itr2->_edge, direction, (DISTRIBUTION_TYPES::enum_t) dist);
#                         map<string, pair<double, double> >::const_iterator f_itr = _statistics_table.find(key);
#                         pair<double, double> stats = f_itr->second;
#                         //cout<<"\t"<<key<<"\t"<<stats.first<<"\n";
#                         if (stats.first>max_stats){
#                             max_stats = stats.first;
#                         }
#                         if (stats.first<min_stats){
#                             min_stats = stats.first;
#                         }
#                     }
#                     min_join_selectivity = min_join_selectivity * min_stats;
#                     max_join_selectivity = max_join_selectivity * max_stats;
#                 }
#             }
#         }
#     }

#     sort(degree_array.begin(), degree_array.end());

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

#         /// Randomly choose variables, ignore types date, integer, name or string...
#         /// You should also ignore variables corresponding to join vertices...
#         vector<string> eligible_list;
#         for (map<string, set<string> >::iterator itr=variable_map.begin(); itr!=variable_map.end(); itr++){
#             string var_type = *(itr->second.begin());
#             if (var_type.compare("date")==0 || var_type.compare("integer")==0 || var_type.compare("name")==0 || var_type.compare("string")==0){
#                 continue;
#             }
#             eligible_list.push_back(itr->first);
#         }
#         random_shuffle(eligible_list.begin(), eligible_list.end());

#         /// Select <const-count> variables...
#         if (eligible_list.size()>const_count){
#             eligible_list.erase(eligible_list.begin()+const_count, eligible_list.end());
#         }

#         bool varValid = true;
#         string mapping = "";
#         for (int vid=0; vid<eligible_list.size(); vid++){
#             string var_name = eligible_list[vid];
#             string var_root = var_name.substr(1);

#             mapping.append("#mapping");
#             mapping.append(" ");
#             mapping.append(var_root); /// Drop preceding '?'...
#             mapping.append(" ");
#             mapping.append(*(variable_map[var_name].begin()));
#             mapping.append(" ");
#             mapping.append("uniform");
#             mapping.append("\n");

#             boost::algorithm::replace_all(query_str, var_name + string(" "), string("%") + var_root + string("%") + string(" "));
#             boost::algorithm::replace_all(query_str, var_name + string("\t"), string("%") + var_root + string("%") + string("\t"));

#             if (!constJoinVertexAllowed){
#                 string placeholder = string("%") + var_root + string("%");
#                 if (query_str.find(placeholder)!=string::npos && query_str.find(placeholder)!=query_str.rfind(placeholder)){
#                     varValid = false;
#                 }
#             }
#         }

#         bool qValid = false;
#         string qTemplate = "";
#         qTemplate.append(mapping);
#         qTemplate.append("SELECT ");
#         for (set<string>::iterator itr=variable_set.begin(); itr!=variable_set.end(); itr++){
#             string var_name = *itr;
#             if (find(eligible_list.begin(), eligible_list.end(), var_name)==eligible_list.end()){
#                 qTemplate.append(var_name);
#                 qTemplate.append(" ");
#                 qValid = true;
#             }
#         }
#         qTemplate.append("WHERE {");
#         qTemplate.append("\n");
#         qTemplate.append(query_str);
#         qTemplate.append("}");
#         qTemplate.append("\n");
#         qTemplate.append("#end");
#         qTemplate.append("\n");

#         if (varValid && qValid){
#             cout << qTemplate;
#             return true;
#         }
#     /*
#     }
#     */

#     return false;
# }
        pass

# void statistics::print_graph() const{
    def print_graph(self):
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
