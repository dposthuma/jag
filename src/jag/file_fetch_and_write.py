"""
@author: maartenk, esther lips

file handling script
"""

import jag.common as common
import os
import shutil
import re
import glob
import time
import logging as log
import sys
import tempfile


class InAndOut(object):
    """
    class to handle reading and writing of files
    """
    def __init__(self):
        self.out = os.path.abspath(os.path.curdir) + "/"
        self.time_extention = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
        self.run_rproject = True
       
    
    def set_outfile(self, out_str):
        """
        set a path prefix in the same way plink --out functions
        """
        if (not out_str == None):
            full_path = os.path.abspath(out_str.strip())                #get full pathname
            full_dir = os.path.dirname(full_path)
            #check dir exists and create if not exists and check if dir is writable
            if (_create_output_dir(full_dir) and _check_writability(full_dir)):
                #if only a directory is set, do not include a zero in out path
                if (os.path.basename(full_path) == ""):
                    self.out = full_path
                else:   
                    self.out = full_path + "."
            else:
                print("Jag has stopped.")
                sys.exit(1)
    
        return True

    
    def save_sumlog(self, assoc_file, clusterresults):
        """
        save sumlog file to disk.
        """
        
        output_number = get_number_of_assoc_file(assoc_file)
        output_number= str(output_number[1:])
        
        filename = self.create_assoc_analysis_filename(output_number)
        results_as_text = clusterresults.format_results()
                
        if os.path.exists(filename):
            _remove_file([filename2 for filename2 in self.get_sumlog_filenames()])

        common.save_text_file(filename, results_as_text)
        log.info("Saved sumlog file as " + filename)


    def copy_assoc_file(self, assoc_file):
        """
            copy assoc file to directory to where program is ran from
    
        """
        #get base extension e.g. qassoc,assoc assoc.linear 
        base_extention = common.plink_results_extension(assoc_file)
        output_number = get_number_of_assoc_file(assoc_file)
        assoc_files = []
        extentions = ["", ".adjusted"]
        
        for ext in extentions:
            if (os.path.exists(assoc_file + ext)):
                dest_assoc_file = self.out + "results" + output_number + "." + base_extention + ext
                
                if os.path.exists(dest_assoc_file):
                    _remove_file(dest_assoc_file)
                shutil.copy(assoc_file, dest_assoc_file)
                os.path.abspath(os.path.curdir)
                log.info("\nResults of PLINK saved as " + dest_assoc_file)
                assoc_files.append(dest_assoc_file)
            
            else:
                pass
            
        return(assoc_files)


    def copy_prune_in_file(self, temp_file):
        """
            copy prun.in file to directory to where program is ran from
    
        """
        prune_file_name = temp_file + ".prune.in"
        
        if (os.path.exists(prune_file_name)):
            dest_file = self.out + "prune.in"
            
            if os.path.exists(dest_file):
                _remove_file(dest_file)
            shutil.copy(prune_file_name, dest_file)
            os.path.abspath(os.path.curdir)
            
        else:
            log.info(prune_file_name + " does not exists")
            
        return(dest_file)


    def create_assoc_analysis_filename(self, number):
        """
        create a full path for a assoc analysis file based on number
        """
        
        return self.out + str(number) + ".sumlog"


    def get_sumlog_filenames(self):
        """
        get association analysis filename from the current directory
        with a format like *.sumlog
        """

        files = glob.glob(self.out + "*.sumlog")
        test = re.compile("(\.[Pp]\d{0,3}){0,1}\.sumlog$")
        files_clean = [x for x in files if test.search(x) != None]

        return files_clean


    def save_text_to_filename(self, filename, text):
        """save random draws """
        filename = self.out + filename
        
        if os.path.exists(filename):
            _remove_file(filename)
            
        common.save_text_file(filename, text)
        
        return(filename)


def get_assoc_file_number(resultfile, assoc_file):
    """
    get the number of the assoc file. If not found return 1
    result file is the base temp file 
    assoc it the whole path of the files
    """
    
    resultfile = re.escape(resultfile)
    match = re.search(resultfile + "(\.P(?P<number>\d{0,3})){0,1}\.(q){0,1}assoc", assoc_file)

    number = match.group("number")

    if (number == None):
        number = 1
        
    return number


def get_number_of_assoc_file(assoc_file):
    """
    get file number from assoc file
    """
    
    match = re.search("(?P<number>\.P\d{0,3}){0,1}\.(q){0,1}assoc*", assoc_file)
    output_number = match.group("number")
    
    if(output_number == None):
        output_number = ".P1"
        
    return output_number


def _remove_file(path):
    """
    removes path and gives message
    input can be a string with path or a list with paths
    """
    if type(path).__name__ == "list":
        for path_item in path:
            os.remove(path_item)

    else:
        os.remove(path)


def _check_writability(dir_path):
    """
    check if a directory is writable by writing a temp file to a file.
    This file is  removed 
    returns true is is writable
    """
    try:
        temp_file = tempfile.TemporaryFile(dir=dir_path)
        temp_file.write('This file can be deleted.This is used by JAG to check\
            writability of a directory')
        temp_file.close()
    except OSError:
        print ("directory not writable")
        return False
    else:
        return True
    
    
def _create_output_dir(full_dir):
    """
    create dir if not already exists
    """
    
    if not os.path.isdir(full_dir):
        try:
            print("Creating the following directory: " + full_dir)
            os.makedirs(full_dir)
            
        except OSError:
            print("Could not create directory: " + full_dir + " ")
            return False
        
    return True

