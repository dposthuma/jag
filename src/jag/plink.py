'''
Created on Dec 16, 2010

@author: maartenk, esther lips

Call the Plink program and return raw results as pathname or as score 
'''

import logging as log
import jag.common as common
import jag.permutationresults as permutationresults
import os
import platform
import re
import subprocess
import sys
import tempfile
import time

class Plink(object):
    """
    handles input and output from plink
    """
    def __init__(self):
        self.bfile = None
        self.pheno_file = None
        self.covar_file = None
        self.plinkcommand = ""
        self.snp_to_group_map = ""
        self.group = ""
        self.arguments = None
        self.switches = " "
        self.phenotypes_present = []
        self.verbose = False
        self.print_command = False
        self.permutations_completed = 0


    def set_plink_arguments(self, arguments):
        """
        set arguments for plink except bfile and output
        """
        if(len(arguments) > 2):
            self.arguments = arguments


    def set_snp_to_group_map(self, snp_to_group_map):
        """
        check snp_to_group_map is accessible for the program
        """
        self.snp_to_group_map = snp_to_group_map


    def get_number_of_pheno_types(self):
        """
        Get the number of phenotypes based on the 
        number of columns in the phenofile
        """
        if (self.pheno_file != None):
            pheno = common.read_pheno_file(self.pheno_file)
            total_columns = len(pheno[1][1].split("\t"))
        else:
            total_columns = 1
        return (total_columns)


    def _print_assoc_state(self, process, nr_of_phenotypes):
        """
        filter output when running a association test in Plink
        """
        output_of_plink = ""
        test = re.compile("\.P(\d{0,4}){0,1}\.(q){0,1}assoc\s")
        for line in iter(process.stdout.readline, ''):
            console_out = line.rstrip()
            match = test.search(console_out)
            
            if (match != None):
                column_number = int(match.groups()[0])
                if (nr_of_phenotypes != 0):
                    if (column_number % nr_of_phenotypes == 0):
                        perm_this_run = column_number / nr_of_phenotypes
                        permutation_completed = str(self.permutations_completed + perm_this_run)
                        log.info("permutation " + permutation_completed + " ready")
                
            output_of_plink += str(console_out) + "\n"

        return output_of_plink


    def run_plink(self):
        """
        run a instance of plink
        """
        outfile = tempfile.mkstemp()

        # select logistic test
        if (self.switches.count("--logistic") >= 1):
            self.arguments = self.arguments.replace("--assoc", "")
            #remove adjust when logistic and covar is set
            if (self.arguments.count("--covar") >= 1):
                self.arguments = self.arguments.replace("--adjust", "")
        
        # select linear test
        if (self.switches.count("--linear") >= 1):
            self.arguments = self.arguments.replace("--assoc", "")
            #remove adjust when linear and covar is set
            if (self.arguments.count("--covar") >= 1):
                self.arguments = self.arguments.replace("--adjust", "")


        command = "plink --noweb " + self.arguments + self.switches + " --out " + outfile[1]
      
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        
        if (self.print_command):
            log.info(self.arguments + self.switches)
   
        output_of_plink = ""
        nr_of_phenotypes = len(self.phenotypes_present) # length is zero
    
        #nice output for assoc files
        if(re.search("--assoc", self.arguments)):
            output_of_plink = self._print_assoc_state(process, nr_of_phenotypes)
        elif(re.search("--linear", self.switches)):
            output_of_plink = self._print_assoc_state(process, nr_of_phenotypes)
        elif(re.search("--logistic", self.switches)):
            output_of_plink = self._print_assoc_state(process, nr_of_phenotypes)
        elif(re.search("--indep-pairwise", self.arguments)):
            output_of_plink = _print_indep_pairwaise_state(process)
        else:
            log.info("no output filtering found for command" + self.arguments)
            for line in iter(process.stdout.readline, ''):
                log.info(line.strip())
                output_of_plink += line
       
        if (self.verbose):
            log.info(output_of_plink)
        #check if plink is finished correctly by checking the output of plink "Analysis finished"
        if not re.search("Analysis finished", output_of_plink):
            log.info("PLINK did not finish analysis. Job is terminated.")
            log.info(output_of_plink)
            sys.exit(1)
        return(outfile[1])



    def point_outputfile_to_phenotype(self, seeds, resultfile, results_extension):
        """
        map the files from the permutation phenotypes to the right permutation 
        """
        resultfiles = {}
        nr_of_phenotypes = len(self.phenotypes_present)

        for phenotype in self.phenotypes_present:
            phenotype=int(phenotype)
            resultfiles[phenotype] = []
            range_until=(len(seeds)*nr_of_phenotypes)
            #the second number in range excludes the range : it does not includes the border
            if (phenotype==nr_of_phenotypes):
                range_until+=1
          
            numbers=range(phenotype,range_until,nr_of_phenotypes)
            resultfiles[phenotype]= [resultfile+".P"+str(number)+"."+results_extension for number in numbers]
            
        
        return resultfiles


    def run_permutation(self, seeds):
        """
        performs permutations  
        """
        permuted_pheno_file = common.create_pheno_permutation_file(self, seeds)
        pheno_arg = " --pheno " + str(permuted_pheno_file + " --all-pheno" + " --extract " + self.group)

        if self.covar_file:
            permuted_covar_file = common.create_covar_permutation_file(self, seeds)
            pheno_arg += " --covar " + permuted_covar_file + " --hide-covar"
            self.set_plink_arguments("--bfile " + self.bfile + pheno_arg)
        else:
            self.set_plink_arguments("--bfile " + self.bfile + " --assoc " + pheno_arg)


        resultfile = self.run_plink() # run plink
        results_extension = common.plink_results_extension(resultfile)
        resultfiles = self.point_outputfile_to_phenotype(seeds, resultfile, results_extension)
        results = extract_permutated_scores(self.snp_to_group_map, resultfiles, seeds)
                        
        os.remove(permuted_pheno_file)  #remove permuted pheno file 
        common.remove_plink_output(resultfile)  # remove all output from plink
     
        return results


    def run_permutation_single_core(self, jobmap):
        """
        run plink on a single cpu
        """
        permutations = []
        
        for seeds in jobmap:
            permutations.append(self.run_permutation(seeds))
            self.permutations_completed += len(seeds)
   
        return permutations


def extract_permutated_scores(snp_to_group_map, resultfiles, seeds):
    """
    extract permuted scores from a assoc files
    """
    permuted_results = {}
    
    for pheno, files in resultfiles.iteritems():
        clusterresults = permutationresults.PermutatedResults()
        
        for resultfile in files:
            filehandle = common.getfile_handle(resultfile)
                        
            for line in filehandle:
                splitted = line.split()
                snp_id = splitted[1]
                p_value = splitted[8]
               
                if(p_value == "NA"):
                    p_value = 1

                if (snp_id in snp_to_group_map):
                    clusterresults.add_snp_to_cluster(p_value, snp_to_group_map[snp_id])

            clusterresults.process_permutation()

        clusterresults.add_seeds(seeds)
        permuted_results[pheno] = clusterresults
        
    return permuted_results



def _print_indep_pairwaise_state(process):
    """
    filter output when running an indep_pairwaise in Plink
    """
    output_of_plink = ""
    test = re.compile("\For chromosome (\d*),")
    for line in iter(process.stdout.readline, ''):
        console_out = line.rstrip()
        match = test.search(console_out)
        if (match != None):
            log.info ("Pruned independent SNPs on chromosome " + match.groups()[0] + " ")

        output_of_plink += str(console_out) + "\n"

    return output_of_plink
