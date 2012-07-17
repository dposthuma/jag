'''
Created on Dec 15, 2010

@author: maartenk, esther lips

Performs permutation over results from AssociationAnalysis (or creates it when not found)
permutation are only done for pheno types where sumlog.P{number}.log is found. 
'''

from association_analysis import AssociationAnalysis
import jag.common as common
import jag.permutationresults as permutationresults
import random
import sys
import logging as log


class Permutation(object):
    '''
    handle the permutations on a group file
    '''

    def __init__(self, plink, group, in_out, no_emp=False):
        self.plink = plink
        self.group = group
        self.files = {}
        self.in_out = in_out
        self.no_emp = no_emp
        self.gene_based = False


    def get_assoc_analysis(self):
        """
        get association analysis file
        """

        if not self.no_emp:
            assoc_analysis = AssociationAnalysis(self.in_out)
            aa_results = assoc_analysis.get_aa_results(self.plink, self.group, self.gene_based)

            self.plink.phenotypes_present = _get_aa_phenotypes_present(aa_results)
        else:
            aa_results = self.plink.get_number_of_pheno_types() * [False]
            self.plink.phenotypes_present = range(1, len(aa_results) + 1)
        return aa_results
    

    def localpermutation(self, nperm, random_seed, gene_group=False):
        '''
        Run permutation on local computer
        '''
        
        self.gene_based = gene_group

        aa_results = self.get_assoc_analysis()
        #create a random seed
        if (random_seed == None):
            random_seed = _create_random_seed(8)

        snp_to_group = common.map_snp_to_geneset(self.group, gene_group)
        self.plink.set_snp_to_group_map(snp_to_group)

        nperm = int(nperm)
        
        permperjob = int(50)

        if self.plink.covar_file:
            permperjob = 1

        jobs = _create_seeds_for_jobs(nperm, permperjob, random_seed)
     
        resultfiles = self.plink.run_permutation_single_core(jobs)                        
        results_merged = merge_results(resultfiles)
               
        for key, perm_score_per_group in results_merged.iteritems():
            perm_out_string = perm_score_per_group.format_permout()
            perm_out_filename = str(random_seed) + ".P" + str(key) + ".perm"
            #save permutation score as table with random seed in filename
            perm_filename = self.in_out.save_text_to_filename(perm_out_filename, perm_out_string)
            log.info("\nSaved permutation results as " + perm_filename)
            key = str(key)

            self.files[key] = {"perm":perm_filename}
                                
            if (self.no_emp == False):
                column_aaresults_mapping = {}
                
                for i in range(len(aa_results)):
                    column_aaresults_mapping[str(aa_results[i].column_number)] = i
                    
                current_aa_results = aa_results[column_aaresults_mapping[key]]
                empp_out_text = perm_score_per_group.format_permutated_results(current_aa_results)

                emp_filename = "P" + key + ".empp"
                emp_filename = self.in_out.save_text_to_filename(emp_filename, empp_out_text)
                log.info("Saved empirical pvalues as " + emp_filename)
                self.files[str(key)] ["empp"] = emp_filename
                
                if self.in_out.run_rproject:
                    import jag.plot_with_r as plot_with_r
                    plotter = plot_with_r.call_r(self.in_out)
                    plotter.draw_dist_plot(self.files, key)


def select_aa_on_pheno(aa_result_list, pheno_number):
    """
    return from a list of aa_results a aa_result with column number "pheno_number"
    """
    for aa_result in aa_result_list:
        if(int(aa_result.column_number) == int(pheno_number)):
            return(aa_result)
    sys.exit("selection of phenotype failed")
    

def merge_list_of_dict_of_list(resultfiles):
    """
    Merge a list of dicts which contains list2, append list2 on bases of dicts keys
    """
    perm_score_per_group = {}
    for key in resultfiles[0].keys():
        perm_score_per_group[key] = []

    for res_dict in resultfiles:
        for key, value in res_dict.iteritems():
            perm_score_per_group[key] = perm_score_per_group[key] + value

    return perm_score_per_group


def merge_results(resultfiles):
    """
    Merge a list of dicts which contains list2, append list2 on bases of dicts keys
    """
    
    perm_per_pheno = {}
    #loop if it is run on multiple cpu
    for result in resultfiles:
        for pheno_nr in result:
            fg_cluster = result[pheno_nr]
            
            if not perm_per_pheno.has_key(pheno_nr):
                perm_per_pheno[pheno_nr] = permutationresults.PermutatedResults()
                
            for groupname in fg_cluster.permutated_scores.keys():
                
                scores = fg_cluster.permutated_scores[groupname]
                perm_per_pheno[pheno_nr].add_results(groupname, scores)
    
            perm_per_pheno[pheno_nr].add_seeds(fg_cluster.seeds)

    return perm_per_pheno


def _create_seeds_for_jobs(nperm, permpercpu, user_set_random_seed):
    """
    create seeds for each permutation job
    """
    random.seed(user_set_random_seed)
    seeds = nperm * [0]
    for i in range(len(seeds)):
        seeds[i] = random.random()
    amountofruns = int(nperm / permpercpu)

    n_perm_non_fair_dist = int(nperm % permpercpu)
    #create first job including non fair permutations

    if (amountofruns == 0):
        amountofruns = 1
    jobs = [0] * amountofruns
    start_of_slice = permpercpu + n_perm_non_fair_dist
    jobs[0] = seeds[:start_of_slice]
    for run in range(1, amountofruns):
        jobs[run] = seeds[start_of_slice:start_of_slice + permpercpu]
        start_of_slice = start_of_slice + permpercpu

    return jobs


def _create_random_seed(length):
    """
    creates a random seed that is a valid filename.
    The length of the random seed can be given as int
    """
    randomized = random.SystemRandom()

    digits = "".join([str(x) for x in range(0, 10)])
    small = "".join([chr(x) for x in range(97, 123)])
    big = "".join([chr(x) for x in range(65, 91)])
    others = "-_"

    all_characters = digits + small + big + others
    return "".join([randomized.choice(all_characters) for x in range(length)])


def _get_aa_phenotypes_present(aa_results):
    phenotypes_present = sorted([int(res.column_number) for res in aa_results])
    return phenotypes_present
