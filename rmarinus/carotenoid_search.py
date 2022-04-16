import numpy, datetime, cobra, multiprocessing, pandas, pickle

#import warnings
#warnings.filterwarnings("ignore")

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

# usage: time python carotenoid_search.py &> carotenoid_search_messages.txt

# 100 x 100 genes = 4,950 pairs. @ M1 takes 21 sec
# 250 x 250 genes = 31,125 pairs. @ M1 takes 2 min
# 578 x 578 genes = 166,753 pairs. @ M1 should take 11 min. It took 12

# 6. parallel computation
if __name__ == '__main__':

    # 0. user-defined variables
    number_of_threads = 8

    # 1. load the model
    model = cobra.io.read_sbml_model("/Users/adrian/scratch/rmarinus/Rmarinus_578_model.xml")
    number_of_genes = len(model.genes)

    # 2. solve the model
    wt_solution = model.optimize()

    # 3. add demand reaction
    model.add_boundary(model.metabolites.get_by_id("CAROT_RMAR_c"), type="demand")

    # 4. growth-coupled design
    reaction_of_interest = 'DM_CAROT_RMAR_c'
    biomass_reaction_label = 'BIOMASS'

    # 5. define tasks
    printt('working with {} genes'.format(number_of_genes))
    tasks = []
    for i in range(len(model.genes)):
        for j in range(len(model.genes)):
            if i < j:
                task = [i, j, reaction_of_interest, biomass_reaction_label, model]
                tasks.append(task)
    printt('working with {} gene pairs'.format(len(tasks)))

    # 6. run in parallel
    with multiprocessing.Pool(number_of_threads) as hydra:
        hydra_output = hydra.map(growth_coupled_analysis, tasks)
    printt('obtained {} results'.format(len(hydra_output)))

    # 7. store results as a dataframe, then pickle it
    df = pandas.DataFrame(hydra_output, columns=['i', 'j', 'KO growth', 'min production', 'max production'])
    df.sort_values(by=['min production'], ascending=False, inplace=True)
    print(df.head())

    printt('store double KO information dataframe as a pickle')
    f = open('doubleKO.pickle','wb')
    pickle.dump(df, f)
    f.close()
