'''
Created on Dec 13, 2010

@author: maartenk, esther lips

Script performs association analysis/self-contained test.

Saves result from plink as results.{.PhenotypeNumber}.assoc 

'''

import logging as log
import glob
import jag.clusterresults as CR
import jag.common as common
import jag.file_fetch_and_write as in_out
import os
import jag.plot_with_r as plot_with_r
import sys



class AssociationAnalysis(object):
    '''
    Class to Association Analysis methods.
    '''

    def __init__(self, inout):
        self.in_and_out = inout
        self.files = {}


    def get_aa_results(self, plinkin,geneset, gene_based=False):
        """
        check if there is already sumlog*.log in current dir and return as cluster 
        results file. If files are not present create the files

        Returns list of cluster results objects
        """
        results_files = self.in_and_out.get_sumlog_filenames()
        aa_results = []
        
        if (len(results_files) > 0):
            for i in range(len(results_files)):
                if(os.path.exists(results_files[i])):
                    file_handle = common.getfile_handle(results_files[i])   #read file
                    aa_result = CR.Clusterresults()
                    aa_result.read_formated_results(file_handle)
                    
                    if (aa_result.geneset_path == geneset):
                        aa_results.append(aa_result)
                    else:
                        log.info("sumlog file does not have identical path to geneset file\n"
                            + "Will run Association analysis again")
                        aa_result = self.run_step1(plinkin, geneset, gene_based)
                        clusterresults = [pheno_dict["clusterresults"] for  pheno_dict in aa_result.values()]
                        aa_results.extend(clusterresults)
                        break
                else:
                    aa_result = self.run_step1(plinkin, geneset, gene_based)
                    clusterresults = [pheno_dict["clusterresults"] for pheno_dict in aa_result.values()]
                    aa_results.extend(clusterresults)
        else:
            clusterresults = [pheno_dict["clusterresults"] for  pheno_dict in self.run_step1(plinkin, geneset, gene_based).values()]
            aa_results = clusterresults

        return aa_results


    def run_step1(self, plink, geneset, gene_group=False, adjusted=False):
        """
        main function of association analysis. this function calls all function 
        needed for aa analysis.
        plink. a plink object. Example to create this kind of object can be found in test/test_jag.py
        group: path to file which contains file mapping 
        """
        #check files are available
        #read geneset as dictionary (key=snp, value=geneset_name)
        if (not geneset == None):
            snptogeneset = common.map_snp_to_geneset(geneset, gene_group)
                    
        pheno_arg = ""
        
        if (not plink.pheno_file == None):
            pheno_arg = " --pheno " + str(plink.pheno_file + " --all-pheno ")
            
        if not plink.covar_file == None:
            pheno_arg += " --covar " + plink.covar_file + " --hide-covar "
            plink.set_plink_arguments("--bfile " + plink.bfile + " " + pheno_arg)
        
        else:
            plink.set_plink_arguments("--bfile " + plink.bfile + " --assoc " + pheno_arg)

        resultfile = plink.run_plink()
        assoc_files = glob.glob(resultfile + "*.*assoc*")
                
        if adjusted:
            assoc_files = [assoc_file for assoc_file in assoc_files if assoc_file.endswith("adjusted")]
                
        else:
            assoc_files = [assoc_file for assoc_file in assoc_files if not assoc_file.endswith("adjusted")]
          
        pheno_nr=1
        for assoc_file in assoc_files: 
            clusterresults = CR.Clusterresults()            #    read values from assoc assoc_file and lookup in snptogroup map:
            clusterresults.pheno_path = plink.pheno_file
            clusterresults.total_columns = plink.get_number_of_pheno_types()
            clusterresults.column_number = str(in_out.get_assoc_file_number(resultfile, assoc_file))
            self.files[clusterresults.column_number] = {}
            clusterresults.bfile_path = plink.bfile
            clusterresults.geneset_path = geneset

            moved_assoc_file = self.in_and_out.copy_assoc_file(assoc_file)
            
            if (not geneset == None):
                clusterresults.map_p_values_from_assoc_file(snptogeneset, assoc_file, adjusted)      # extract need p values from assoc file
                self.files[clusterresults.column_number]["raw_p_files"] = clusterresults.write_all_pvalues()
                clusterresults.createcleanresults()
                self.files[clusterresults.column_number]["clusterresults"] = clusterresults
                self.in_and_out.save_sumlog(assoc_file, clusterresults) #'generate .sumlog'
            
            else:
                log.info ("No geneset file is set: no sumlog file is made.")
                sys.exit("JAG is terminated.")
            
            # code for QQ plots    
            self.files[clusterresults.column_number]["assoc_files"] = moved_assoc_file              #saving full file path for creating plot later on
            
            assoc_pheno = [values["assoc_files"] for values in self.files.itervalues()]           
            assoc_files = [ assoc_file for assoc_file in plot_with_r.from_iterable(assoc_pheno)]
            
            raw_p_pheno = [values["raw_p_files"] for values in self.files.itervalues()] 
            set_files = [ assoc_file for assoc_file in plot_with_r.from_iterable(raw_p_pheno)]
            
            if self.in_and_out.run_rproject:
                plotter = plot_with_r.call_r(self.in_and_out)
                plotter.draw_qq_plots(assoc_files, pheno_nr, adjusted)
                plotter.draw_qq_plots_for_sets(raw_p_pheno, pheno_nr)
                
            pheno_nr=pheno_nr+1
                       
        common.remove_plink_output(resultfile)          # remove all output from plink
        
        return self.files
