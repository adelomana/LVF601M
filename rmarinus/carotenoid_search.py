import numpy, pandas, datetime
import cobra, cobra.test

import multiprocessing, multiprocessing.pool
from multiprocessing import Process, Queue

def growth_coupled_analysis(task):

    """
    This function performs the growth-coupled production.
    It takes as input a list as [first_gene_pair_index, second_gene_pair_index, reaction_of_interest, biomass_reaction_label]
    It gives as output a list as [first_gene_pair_index, second_gene_pair_index, growth, min_production, max_production]
    """

    i = task[0]
    j = task[1]
    reaction_of_interest = task[2]
    biomass_reaction_label = task[3]
    model = task[4]

    with model as model:

        # KO
        model.genes[i].knock_out()
        model.genes[j].knock_out()
        solution = model.optimize()
        if solution.status == 'optimal':
            ko_growth = solution.objective_value

            # growth-coupled production
            model.objective = reaction_of_interest
            model.reactions.get_by_id(biomass_reaction_label).lower_bound = ko_growth
            max_production = model.optimize(objective_sense='maximize').objective_value
            min_production = model.optimize(objective_sense='minimize').objective_value

            result = [i, j, ko_growth, min_production, max_production]
        else:
            result = [i, j, 0, 0, 0]

    return result

def printt(message):

    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(message)))

    return None

###
### MAIN
###

# 1. load model
model = cobra.io.read_sbml_model("/Users/adrian/scratch/rmarinus/Rmarinus_578_model.xml")
print(model.summary())

# 2. get info about model
wt_solution = model.optimize()
print(wt_solution.objective_value)

number_of_genes = len(model.genes)
print(number_of_genes)

# 3. create demand reaction
model.add_boundary(model.metabolites.get_by_id("CAROT_RMAR_c"), type="demand")
print("demands", model.demands)
for reaction in model.demands:
    print(reaction)

# 4. growth-coupled design
reaction_of_interest = 'DM_CAROT_RMAR_c'
biomass_reaction_label = 'BIOMASS'

# 4.1. run in a paralel environment
number_of_threads = 2

printt('working with {} genes'.format(number_of_genes))

tasks = []
for i in range(len(model.genes[:5])):
    for j in range(len(model.genes[:5])):
        if i < j:
            task = [i, j, reaction_of_interest, biomass_reaction_label, model]
            tasks.append(task)
printt('working with {} gene pairs'.format(len(tasks)))

printt('entering a parallel world of {} threads'.format(number_of_threads))
hydra = multiprocessing.pool.Pool(number_of_threads)
hydra_output = hydra.map(growth_coupled_analysis, tasks)
hydra.close()
printt('completed {} tasks'.format(len(hydra_output)))
