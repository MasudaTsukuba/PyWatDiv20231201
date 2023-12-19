"""Generate SPARQL queries
query_generator.py
Generate SPARQL queries from the source RDF file.
Queries give at least one result for the source RDF file.
T. Masuda, 2023/12/18
"""

import random
import sys
from rdflib import Graph, Variable


class Vertex:
    """
    Class for handling a vertex.
    Vertex is either a subject or an object.

    Attributes:
        left_vertices (list[Vertex]): list of vertices that become subjects of triples while the current vertex becomes an object.
        right_vertices (list[Vertex]): list of vertices that become objects of triples while the current vertex becomes a subject.
        uri (str): uri string of this vertex with < and >
        predicate (str): predicate string from parent (either subject or object) to this vertex
        is_variable (bool): whether this vertex becomes a variable or not
    """

    var_number: int = 0  # integer for generating a name of a variable, such as ?v0.

    def __init__(self, graph: Graph):
        """
        Constructor of Vertex class.
        """
        self.left_vertices: list[Vertex] = []
        self.right_vertices: list[Vertex] = []
        self.uri: str = ''  # uri of this vertex
        self.predicate: str = ''  # predicate from parent to this vertex
        self.is_variable: bool = False  # whether this vertex becomes a variable or not
        self.graph: Graph = graph  # RDF graph
        Vertex.var_number = 0  # reset the variable number

    def expand_vertex(self) -> bool:
        """Expand vertices starting from the root vertex.
        Expansion info is stored in left_vertices and right_vertices of each vertex.

        Returns:
            bool: expansion is successful
        """
        def find_vertex(side='right'):
            """Find a vertex connected to the current vertex.

            Args:
                side (str): which side, left or right, the vertex is searched

            Returns:

            """
            # q = ""
            if side == 'right':
                if self.uri.find('http://') < 0:  # uri is a literal
                    return False, []  # return with failure
                q = f"""
                    SELECT DISTINCT ?p ?v WHERE {{ {self.uri} ?p ?v .}}
                    """
            else:  # right side (=object) can be a literal
                q = f"""
                    SELECT DISTINCT ?p ?v WHERE {{ ?v ?p {self.uri} .}}
                    """
            # print(q)  # debug
            print('.', end='')  # debug for showing the progress
            results = None  # no results are assumed
            try:
                results = self.graph.query(q)  # execute the query
            except Exception as e:
                pass  # some trouble in the query
            if len(results.bindings) == 0:  # no results are found
                return False, []  # return with failure
            list_of_vertices = []
            for result in results.bindings:
                vertex = str(result[Variable('v')])
                edge = str(result[Variable('p')])
                if edge.find('http://') >= 0:  # predicate must be an uri
                    list_of_vertices.append((vertex, edge))  # save as a tuple
            random.shuffle(list_of_vertices)  # random shuffle
            return True, list_of_vertices

        def expand_left():
            """Expand the left side.
            Root vertex becomes an object of a triple.

            Returns:

            """
            if self.left_vertices:
                success, list_of_vertices = find_vertex(side='left')
                if success:
                    if len(list_of_vertices) < len(self.left_vertices):
                        return False  # not enough vertices found
                    success = True  # assume success
                    random.shuffle(list_of_vertices)
                    for index, left_vertex in enumerate(self.left_vertices):
                        if list_of_vertices[index][0].find('http://') >= 0:
                            left_vertex.uri = f'<{list_of_vertices[index][0]}>'
                        else:
                            left_vertex.uri = f'"{list_of_vertices[index][0]}"'
                            # left_vertex.uri = f'{list_of_vertices[index][0]}'
                        # left_vertex.predicate = f'<{list_of_vertices[index][1]}>'
                        left_vertex.predicate = f'{list_of_vertices[index][1]}'
                        success_left_vertex = left_vertex.expand_vertex()
                        if not success_left_vertex:
                            success = False
                            break
                return success
            else:
                return True  # no further left side vertices

        def expand_right():
            """Expand the right side.
            Root vertex becomes a subject of a triple.

            Returns:

            """
            if self.right_vertices:
                success, list_of_vertices = find_vertex(side='right')
                if success:
                    if len(list_of_vertices) < len(self.right_vertices):
                        return False
                    success = True
                    random.shuffle(list_of_vertices)
                    for index, right_vertex in enumerate(self.right_vertices):
                        if list_of_vertices[index][0].find('http://') >= 0:
                            right_vertex.uri = f'<{list_of_vertices[index][0]}>'
                        else:
                            right_vertex.uri = f'"{list_of_vertices[index][0]}"'
                            # right_vertex.uri = f'{list_of_vertices[index][0]}'
                        # right_vertex.predicate = f'<{list_of_vertices[index][1]}>'
                        right_vertex.predicate = f'{list_of_vertices[index][1]}'
                        success_right_vertex = right_vertex.expand_vertex()
                        if not success_right_vertex:
                            success = False
                            break
                return success
            else:
                return True  # no further right side vertices

        success_left = expand_left()  # expand the left side
        success_right = expand_right()  # expand the right side
        return success_left and success_right  # both sides must be succeeded

    def generate_query(self):
        """Generate a SPARQL query starting from this vertex.

        Returns:
            str: query string
        """
        query_string = 'SELECT $var_list WHERE { \n$triple_list }'
        list_of_all_variables = []
        list_of_all_triples = []
        subject = self.uri
        if self.is_variable:
            subject = f'?v{str(Vertex.var_number)}'
            Vertex.var_number += 1
            list_of_all_variables.append(subject)
        for left_vertex in self.left_vertices:
            list_of_variables, list_of_triples = left_vertex.generate_query_parts(subject, side='left')
            list_of_all_variables += list_of_variables
            list_of_all_triples += list_of_triples

        for right_vertex in self.right_vertices:
            list_of_variables, list_of_triples = right_vertex.generate_query_parts(subject, side='right')
            list_of_all_variables += list_of_variables
            list_of_all_triples += list_of_triples

        var_list = ' '.join(list_of_all_variables)
        triple_list = '\n'.join(list_of_all_triples)
        query_string = query_string.replace('$var_list', var_list).replace('$triple_list', triple_list)
        return query_string

    def generate_query_parts(self, vertex, side='right'):
        """

        Args:
            vertex:
            side:

        Returns:

        """
        list_of_all_variables = []
        list_of_all_triples = []
        if side == 'right':
            # object_ = f'"{self.uri}"'
            object_ = f'{self.uri}'
            if self.uri.find('http://') >= 0:
                # object_ = f'<{self.uri}>'
                object_ = f'{self.uri}'

            if self.is_variable:
                object_ = f'?v{Vertex.var_number}'
                Vertex.var_number += 1
                list_of_all_variables.append(object_)
            triple = f'{vertex} <{self.predicate}> {object_} . '
            list_of_all_triples.append(triple)
            for left_vertex in self.left_vertices:
                list_of_variables, list_of_triples = left_vertex.generate_query_parts(object_, side='left')
                list_of_all_variables += list_of_variables
                list_of_all_triples += list_of_triples
            for right_vertex in self.right_vertices:
                list_of_variables, list_of_triples = right_vertex.generate_query_parts(object_, side='right')
                list_of_all_variables += list_of_variables
                list_of_all_triples += list_of_triples

        else:
            subject = f'{self.uri}'
            if self.uri.find('http://') >= 0:
                # subject = f'<{self.uri}>'
                subject = f'{self.uri}'

            if self.is_variable:
                subject = f'?v{Vertex.var_number}'
                Vertex.var_number += 1
                list_of_all_variables.append(subject)
            triple = f'{subject} <{self.predicate}> {vertex} . '
            list_of_all_triples.append(triple)
            for left_vertex in self.left_vertices:
                list_of_variables, list_of_triples = left_vertex.generate_query_parts(subject, side='left')
                list_of_all_variables += list_of_variables
                list_of_all_triples += list_of_triples
            for right_vertex in self.right_vertices:
                list_of_variables, list_of_triples = right_vertex.generate_query_parts(subject, side='right')
                list_of_all_variables += list_of_variables
                list_of_all_triples += list_of_triples

        return list_of_all_variables, list_of_all_triples


def generate_pattern(graph: Graph):
    """Generate query branch pattern.

    Args:
        graph (Graph): RDF graph containing triples in *.nt file.

    Returns:
        Vertex: root vertex

    """
    def generate_linear():
        """

        Returns:

        """
        vertex3 = Vertex(graph)
        vertex3.is_variable = True

        vertex2 = Vertex(graph)
        vertex2.right_vertices.append(vertex3)
        vertex2.is_variable = True

        vertex1 = Vertex(graph)
        vertex1.right_vertices.append(vertex2)
        vertex1.is_variable = False

        return vertex1

    def generate_linear_both_sides(num_left: int = 1, num_right: int = 1):
        """

        Args:
            num_left:
            num_right:

        Returns:

        """
        root_vertex = Vertex(graph)
        root_vertex.is_variable = True

        vertex_right = root_vertex
        for i in range(num_left):
            new_vertex = Vertex(graph)
            vertex_right.left_vertices.append(new_vertex)
            vertex_right = new_vertex

        vertex_left = root_vertex
        for i in range(num_right):
            new_vertex = Vertex(graph)
            vertex_left.right_vertices.append(new_vertex)
            vertex_left = new_vertex

        return root_vertex

    def generate_star(num_left: int = 0, num_right: int = 2):
        """Generate star-like query pattern.

        Args:
            num_left (int): number of edges on the left side
            num_right (int): number of edges on the right side

        Returns:
            Vertex: root vertex

        """
        root_vertex = Vertex(graph)
        root_vertex.is_variable = True

        vertex_right = root_vertex
        for i in range(num_left):
            new_vertex = Vertex(graph)
            vertex_right.left_vertices.append(new_vertex)

        vertex_left = root_vertex
        for i in range(num_right):
            new_vertex = Vertex(graph)
            vertex_left.right_vertices.append(new_vertex)

        return root_vertex

    def generate_tree(level1: int = 2, level2: int = 2):
        """Generate tree shaped query pattern.

        Returns:
            Vertex: root vertex

        """
        root_vertex = Vertex(graph)
        root_vertex.is_variable = True
        for i in range(level1):
            new_vertex = Vertex(graph)
            root_vertex.right_vertices.append(new_vertex)
            for j in range(level2):
                new_vertex2 = Vertex(graph)
                new_vertex.right_vertices.append(new_vertex2)
        return root_vertex

    vertex = generate_linear()
    vertex = generate_linear_both_sides(1, 1)  # vertex->start->vertex
    vertex = generate_linear_both_sides(0, 2)  # start->vertex->vertex
    vertex = generate_star(0, 3)  # star shaped query, start->vertex1, start->vertex2
    # vertex = generate_star(3, 0)  # star shaped query, vertex1->start, vertex2->start
    # vertex = generate_tree(2, 0)  # tree shaped query
    return vertex


def list_predicates(graph: Graph):
    """

    :return:
    """

    q = """
    SELECT DISTINCT ?p WHERE { ?s ?p ?o .}
    """

    results = graph.query(q)
    list_of_predicates = []

    for binding in results.bindings:
        list_of_predicates.append(str(binding[Variable('p')]))

    return list_of_predicates


def find_start_point(graph: Graph):
    """

    :return:
    """
    # find a start point (vertex)
    q = """
    SELECT DISTINCT ?s WHERE { ?s ?p ?o .}
    """
    results = graph.query(q)
    list_of_subject = []
    for binding in results.bindings:
        result = str(binding[Variable('s')])
        if result.find('http://') >= 0:
            list_of_subject.append(f'<{result}>')
    if len(list_of_subject) == 0:
        print('No valid subject found')
        sys.exit(-1)
    random.shuffle(list_of_subject)
    start_point = list_of_subject[0]
    return start_point


def execute_query(graph, query):
    """

    :param graph:
    :param query:
    :return:
    """
    results = graph.query(query)
    return results


def main():
    g = Graph()
    g.parse('../SWDF_small20231215.nt')

    # list_of_predicates = list_predicates(g)

    root = generate_pattern(g)  # generate query pattern

    with open('../generated_query.txt', 'a') as output_file:
        repeat = 5  # times of query generation
        for i in range(repeat):
            failed = True
            while failed:
                start_point = find_start_point(g)  # find a vertex that becomes a start point of creating SPARQL query
                root.uri = start_point  # set the actual uri as the start point

                success = root.expand_vertex()  # find an actual sub graph satisfying the specified branch pattern
                if success:
                    print()
                    query_string = root.generate_query()  # generate a SPARQL query string from the expanded vertices pattern
                    print(query_string)
                    output_file.writelines(query_string)
                    output_file.write('\n=============================================================================================================================================================\n')
                    failed = False
                    results = execute_query(g, query_string)  # execute the generated SPARQL query
                    print('Number of query results: ', len(results.bindings))  # and check how may results are obtained
                    Vertex.var_number = 0  # reset the number for ?v0, etc
                    pass  # end of if
                pass  # end of while
            pass  # end of for
    pass  # end of def main():


if __name__ == '__main__':
    main()
