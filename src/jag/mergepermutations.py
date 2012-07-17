'''
Created on Feb 3, 2011

@author: maartenk, esther lips
'''

from jag.clusterresults import Clusterresults
from jag.permutationresults import PermutatedResults
import jag.common as common
import glob
import os
import re
import sys
import logging as log


class MergePermutation:
    """
    Handles the results of all permutations
    """
    def __init__(self, inout):
        self.inout = inout
        self.current_path = os.path.abspath(os.path.curdir) + "/"
        self.files = {}


    def _get_ordered_files(self):
        """
        get all permuted files and return it in a dict with as key phenotype and a list of filenames
        """
        
        files = glob.glob(self.inout.out + "*.perm")   #get perm *.perm
        pattern = re.compile(self.inout.out + ".{1,8}\.P(?P<pheno>.*)\.perm")
        files_clean = [x for x in files if pattern.search(x) != None]
            
        if (len(files) < 1):
            log.info("Could not find any permutations named as " + self.inout.out + "*.perm") 
            log.info(common.get_terminated_time())
            sys.exit()
            
        ordered_results = _get_ordered_files_from_list(files_clean)
        return (ordered_results)


    def _get_ordered_results(self):
        """
        returns dictionary of results by phenotype
        """
        ordered_results = self._get_ordered_files()
        
        for pheno_name in ordered_results.keys():  #read the files
            permutated = PermutatedResults()
            permutated.read_permout(ordered_results[pheno_name])
            ordered_results[pheno_name] = permutated
            
        return ordered_results


    def _load_sumlog_files(self):
        """
        read sumlog files
        """
        
        sumlogfiles = glob.glob(self.inout.out + "*.sumlog")
        test = re.compile("\.P(?P<pheno>.*)\.sumlog$")
        files_clean = [x for x in sumlogfiles if test.search(x) != None]
        aa_result = {}
        
        for file_clean in files_clean:
            aa_result[test.search(file_clean).group("pheno")] = file_clean

        return aa_result


    def mergeresults(self):
        """
        function to merge multiple permutation files into one file
        """
        #load sumlog files to create sumlog files
        aa_result = self._load_sumlog_files()
        ordered_results = self._get_ordered_results()
        perm_files = self._get_ordered_files()
        keys = ordered_results.keys()
        keys.sort(key=common.alpha_sort) 
             
        pheno_nr = 0
        
        #merge each phenotype
        for key in keys:
            pheno_nr=pheno_nr+1
            log.info("\nMerging " + (str(len(perm_files[key]))) + " permutation files for phenotype " + key + "...")
            perm_out_string = ordered_results[key].format_permout() # concatenated results
            perm_filename = "merged.P" + key + ".perm"
            perm_out_filename = self.inout.save_text_to_filename(perm_filename, perm_out_string)
            log.info("Saved merged permutations as " + perm_out_filename)
            #save empirical P file if sumlog files is found
            if (aa_result.has_key(key)):

                aa_object = Clusterresults()
                aa_object.read_formated_results(common.getfile_handle(aa_result[key]))

                empp_out_as_text = ordered_results[key].format_permutated_results(aa_object)
                emp_filename = "merged.P" + key + ".empp"
                empp_filename = self.inout.save_text_to_filename(emp_filename, empp_out_as_text)
                log.info("Saved empirical pvalues as " + empp_filename)
                #call R for distribution plot
                self.files[key] = {"perm":perm_out_filename, "empp":empp_filename}
                
                if self.inout.run_rproject:
                    import jag.plot_with_r as plot_with_r
                    plotter = plot_with_r.call_r(self.inout)
                    plotter.draw_dist_plot(self.files, key)
        
            else:
                log.info("\nWarning: Could not find sumlog file " + self.inout.out + ".P" + key + ".sumlog")
                log.info(common.get_terminated_time())
                sys.exit()
                
                
def _get_ordered_files_from_list(files):
    """
    files:a list of filenames
    return:  a dict with as key phenotype and a list of filenames
    """
    # Skip 'merged' files, but accept manually entered seeds of any length
    test_merged = re.compile("\.merged\.P(?P<pheno>.*)\.perm$")
    files_clean = [x for x in files if test_merged.search(x) is  None]
    
    test = re.compile("\..{1,8}\.P(?P<pheno>.*)\.perm$")
    files_clean = [x for x in files_clean if test.search(x) != None]
    phenotypes = set([test.search(x).group("pheno") for x in files_clean])
    
    ordered_results = {}
    
    for phenotype in phenotypes:
        test_pheno = re.compile("\..{1,8}\.P" + str(phenotype) + "\.perm$")
        pheno_type_files = [x for x in files_clean if test_pheno.search(x) != None]
        ordered_results[phenotype] = pheno_type_files
        
    return ordered_results

