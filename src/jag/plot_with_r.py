'''
Created on Apr 12, 2011

@author: maartenk, esther lips
'''
import os
import subprocess
import tempfile
import logging as log

class call_r(object):
    '''
    Class to group methods to call R.
    '''

    def __init__(self, inout):
        '''Initialize Object to call R, with argument inout to 
        direct produced files to correct output folder.'''
        self.inout = inout        #IO handle

        #Find R binary using 'which'
        self.r_binary = None
        process = subprocess.Popen("which R", shell=True, stdout=subprocess.PIPE)
        exit_status_which = process.wait()
        if exit_status_which == 0:
            self.r_binary = process.communicate()[0].strip()
        else:
            log.info ("R not found" + str(exit_status_which))


    def execute(self, r_script, var_dict, outfile=None):
        '''
        Executing function
        '''
        #Skip if we could not find an R binary
        if not self.r_binary:
            return

        args_text = ""
        if not outfile:
            outfile = tempfile.mkstemp()[1]

        for key, val in var_dict.iteritems():
            args_text += str(key) + "=\"" + str(val) + "\" "
        command = self.r_binary + " CMD BATCH --no-save --no-restore '--args "\
         + args_text + "' " + r_script + " " + outfile

        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        exit_status = process.wait()
        if exit_status == 0:
            pass
        else:
            stdout, stderr = process.communicate()
            log.info("Note: R could not be executed correctly.\n")
            pass


    def draw_dist_plot(self, merged_files, pheno_key):
        '''
        Draw distribution plots
        '''
        files = merged_files[pheno_key]
        log.info("Saved distribution plot as " + self.inout.out + "distribution_sumlogs.P" + pheno_key + ".pdf")
        r_variables = {"emp_file":files["empp"], "out_file":files["perm"], \
                            "prefix" :self.inout.out, "pheno":pheno_key}
        self.execute(_get_package_path() + "distribution_plot.R", r_variables)

   
    def draw_qq_plots(self, assoc_files, pheno_nr, adjusted):
        """Draw Quantile-Quantile plots from assoc files."""
            
        assoc_file = assoc_files[pheno_nr-1]
        image_filename = self.inout.out + "assoc_qq_plot.P" +str(pheno_nr)+ ".pdf"
        log.info("Saved QQ plot from association analysis as " + image_filename)
        r_variables = {"filename":image_filename, "assocfile":assoc_file, \
                           "header":os.path.basename(assoc_file), \
                            "r_path":_get_package_path(), "adjusted":adjusted}
        self.execute(_get_package_path() + "generate_single_qq.R", r_variables)


    def draw_qq_plots_for_sets(self, set_files, pheno_nr):
        """Draw quantile-quantile plots for P values of gene sets."""
        
        set_file = set_files[pheno_nr-1]
        image_filename = self.inout.out + "qq-plot_all_sets.P" + str(pheno_nr) + ".pdf"
        log.info("Saved QQ plot(s) from gene set(s) as " + image_filename)
        
        r_variables = {"p_values_file" : set_file , \
                            "image_filename": image_filename, \
                            "r_path":_get_package_path()}
        self.execute(_get_package_path() + "multiple_qq_plot_for_sets.R", r_variables)


def convert_list_to_arrayray(list_to_convert):
    """
    Convert a python list to a R array (this array is an string in python)
    """
    r_array_text = str(list_to_convert).replace("[", "c(").replace(']', ')')
    return r_array_text


def from_iterable(iterables):
    """
    this function is only provided in python 2.6 and higher
    in the itertools.chain.from_iterable module. 
    This is a drop in replacement for this function
    """
    for iterabl in iterables:
        for element in iterabl:
            yield element
            


def _get_package_path():
    """
    get the absolute path of this file
    """
    return os.path.dirname(os.path.abspath(__file__)) + "/"   
