"""
Test SWDF_small20231215.nt data file
T. Masuda, 2023/12/15
"""

from rdflib import Graph

g = Graph()

g.parse('../SWDF_small20231215.nt')

q = """
SELECT DISTINCT ?p WHERE { ?s ?p ?o .}
"""

results = g.query(q)
pass
