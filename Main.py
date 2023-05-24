#
#  Copyright (c) 2023 - The SUOG Project.
#  @Author: Mirna El Ghosh
#

from core.SimilarityComputationMatrix import SimilarityComputation


def similarity_computation(src, input, measure):

    sc = SimilarityComputation(src, input, measure)

    sc.similarity_images()

    # sc.similarity_phenotypes()


if __name__ == "__main__":

    ontology = 'ontology/suog_ontology_v3.66i1.owl'
    input_file = 'csv/images.csv'
    
    # input_file = 'csv/phenotypes/suog.csv'

    # ontology = 'ontology/hp.owl'
    # input_file = 'csv/phenotypes/hpo.csv'

    similarity_computation(ontology, input_file, "sim_suog") # sim_ic, sim_resnik

