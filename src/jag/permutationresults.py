'''
Created on Dec 14, 2010

@author: maartenk, esther lips
'''

from __future__ import division

import jag.common as common
import sys
import re
import logging as log

class PermutatedResults:
    """
    Handles the results of all permutations
    """
    def __init__(self):
        self.raw_p_value = {}
        self.permutated_scores = {}
        self.seeds = []
        

    def add_seeds(self, seeds):
        """add seeds to permuted results to make reproducing possible"""
        if not isinstance(seeds, list):
            raise AssertionError
        self.seeds.extend(seeds)


    def add_results(self, groupname, scores):
        """
        add score to list of permuted scores . Values are stored per group name
        """
        if not self.permutated_scores.has_key(groupname):
            self.permutated_scores[groupname] = []
        #add array
        self.permutated_scores[groupname].extend(scores)


    def add_snp_to_cluster(self, p_value, gene_and_cluster_list):
        """
        add a P value found in PLINK output to the cluster in gene_and_cluster_list
        """
        for cluster in gene_and_cluster_list:
            if  not cluster["c"] in self.raw_p_value:
                self.raw_p_value[cluster["c"]] = []

            self.raw_p_value[cluster["c"]].append(float(p_value))


    def process_permutation(self):
        """
        calculate the negsumlog10 values from the raw_p_values and add to right cluster in permutated_scores
        remove raw_p_values
        """
        for cluster in self.raw_p_value.iterkeys():
            
            if not cluster in self.permutated_scores:
                self.permutated_scores[cluster] = []
                
            nsumlog10 = common.negsumlog10(self.raw_p_value[cluster])
            self.permutated_scores[cluster].append(nsumlog10)

        #clean the raw p values for next permutation
        self.raw_p_value = {}
        return self.permutated_scores


    def read_permout(self, files):
        """
        read multiple permout files and make one permutationresults file
        """
        for permfile in files:
            file_handle = common.getfile_handle(permfile)
            text_as_list = file_handle.readlines()#pylint: disable=E1103
            header = text_as_list[0].strip().split("\t")
            #check header last column is seed
            if(not header[len(header) - 1] == "seed"):
                sys.exit("last column of " + permfile + " does not have the \"seed\"")
            #create empty list of list to store columns
            results = [  [] for i in range(len(header))]

            for line in range(1, len(text_as_list)):
                splitted_line = text_as_list[line].strip().split("\t")
                for i in range(len(splitted_line)):
                    results[i].append(float(splitted_line[i]))
            #put the results list of list in permutation results
            for i in range(len(header) - 1):
                self.add_results(header[i], results[i])
            self.add_seeds(results[len(header) - 1])

               
    def format_permout(self):
        """format the permutations score in a table like manner with as column headers
        groupname. 
        """

        groupnames = self.permutated_scores.keys()
        groupnames.sort(key=common.alpha_sort)
        length_permutations = [len(self.permutated_scores[g]) for g in groupnames]
       
        diffrence_in_length = sum([abs(l - length_permutations[0])for l in length_permutations])
        formated_results = '\t'.join(groupnames) + "\tseed\n"
        
        if(diffrence_in_length == 0):
            for i in range(length_permutations[0]):
                formated_results += '\t'.join([str(self.permutated_scores[p][i])for p in groupnames])
                formated_results += "\t" + str(self.seeds[i]) + "\n"
        else:
            log.info("\nWarning: gene-sets are not equal over all permuted files. Analysis will be terminated.")
            log.info(common.get_terminated_time())
            exit()
        return(formated_results)


    def format_permutated_results(self, aa_result):
        """
        return nice formated results of permutation with statistics like variance and ngenes etc
        """
        keys = self.permutated_scores.keys()
        keys.sort(key=common.alpha_sort) 
               
        message = "Geneset\tsumlogReal\tnperm\temp_p\tnSNP\tnGenes"
        message += "\tvar(Perms)\tmean(Perms)\tnEff\n"
        for groupname in keys:
            sum_neg_log = aa_result.get_sum_neglog_10_score(groupname)
            n_snp = aa_result.get_n_snps(groupname)
            group_score = self.permutated_scores[groupname]
            message += groupname + "\t"
            message += str(sum_neg_log) + "\t"
            message += str(len(group_score)) + "\t"
            message += str(calculate_emperical_p_value(sum_neg_log, group_score)) + "\t"
            message += str(n_snp) + "\t"
            message += str(aa_result.get_n_genes_of_group(groupname)) + "\t"
            message += str(common.variance(group_score)) + "\t"
            message += str(common.average(group_score)) + "\t"
            message += str(calculate_effective_snps(n_snp, group_score)) + "\n"
        
        return message



def calculate_emperical_p_value(score, permutationscores):
    """
    calculate empirical pvalue - self contained test
    """

    score = float(score)
    permutationscores = [float(x) for x in  permutationscores]
    amount_perm_bigger_score = len([x for x in  permutationscores if x > score])
    
    return (amount_perm_bigger_score / len(permutationscores))



def calc_emp_emp(score, permutationscores):
    """
    calculate empirical pvalue - competetive test 
    """

    score = float(score)
    permutationscores = [float(x) for x in  permutationscores]
    amount_perm_bigger_score = len([x for x in  permutationscores if x < score])
    
    return (amount_perm_bigger_score / len(permutationscores))


def calculate_effective_snps(n_snps, perm):
    """
    in R:
    nEff=trunc(nSNPs*nSNPs*0.189/var(Perms))
    """
    effective_snps = -1
    variance = common.variance(perm)
    if (variance == "NA"):
        effective_snps = "NA"
    elif(variance == 0.0):
        effective_snps = "NA"
    else:
        effective_snps = (int(float(n_snps) * float(n_snps) * 0.189 / variance))
    return(effective_snps)


def read_permutated_results(permfile):
    """
    read the permutated results of a a permutation file and return the 
    header and file as list of lists
    """
    file_handle = common.getfile_handle(permfile)
    text_as_list = file_handle.readline()
    header = text_as_list.strip().split("\t")
    results = [[] for i in range(len(header))]
    for line in file_handle:
        splitted_line = line.strip().split("\t")
        results[0].append(splitted_line.pop(0))
        for i , value in enumerate(splitted_line, 1):
            results[i].append(float(value))

    return (header, results)



    