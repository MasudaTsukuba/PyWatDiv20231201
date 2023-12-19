"""
main_swdf.py
WatDiv applied to SWDF.nt data
T. Masuda, 2023/12/15
"""

from src.model import Model, TripleST, Statistics, QueryTemplateMT


def main(argv: list[str]):
    argc: int = len(argv)  # number of arguments. argv[0] holds the program name.
    model_filename: str = argv[2]  # file name of the model
    cur_model: Model = Model(model_filename)  # read the model

    if argc == 4 and argv[1] == '-d':
        scale_factor: int = int(argv[3])
        cur_model.generate(scale_factor)
        cur_model.save('../data_swdf/saved.txt')  # save randomly generated setting
        return 0
        pass

    if argc == 6 and argv[1] == '-q':
        cur_model.load('../data_swdf/saved.txt')
        query_count = int(str(argv[(argc - 2)]))  # number of queries to be generated
        recurrence_factor = int(str(argv[(argc - 1)]))
        workload = []
        for template_id in range(3, argc - 2):  # file names for query templates
            query_filename = argv[template_id]
            q_template = QueryTemplateMT(mdl=cur_model)  # create an instance of QueryTemplateMT class
            q_template.parse(query_filename)  # read the query template file
            workload = q_template.instantiate(query_count, recurrence_factor)  # receive the generated queries
        for query in workload:
            print(query)
        return 0
        pass

    if argc == 7 and argv[1] == '-s':  # '-s' generates query templates
        cur_model.load('../data_swdf/saved.txt')
        triple_array: list[TripleST] = TripleST.parse_file(argv[3])
        max_q_size: int = int(argv[4])
        q_count: int = int(argv[5])
        const_count: int = int(argv[6])
        Statistics(cur_model, triple_array, max_q_size, q_count, const_count, False, False)
        return 0
    pass


if __name__ == '__main__':
    # argv = ['watdiv', '-d', '../data_swdf/model/swdf-model.txt', '1']  # option -d
    # argv = ['watdiv', '-q', '../data_swdf/model/swdf-model.txt', '../data_swdf/testsuite/test_template.txt', '4', '1']  # option -q
    argv = ['watdiv', '-s', '../data_swdf/model/swdf-model.txt', '../../SWDF_small20231215.nt', '1', '1', '1']  # option -s
    main(argv)

