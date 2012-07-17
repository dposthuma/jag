'''
Created on Dec 13, 2010

@author: maartenk, esther lips

stores all functions which are called from multiple classes
'''

import logging as log
import array
import glob
import math
import os.path
import random
import sys
import tempfile
import re
import time

    

def plink_results_extension(resultfile):
    """
    """
    files_create_by_plink = glob.glob(resultfile + "*")
    return _plink_results_extension(files_create_by_plink)


def _plink_results_extension(list_of_files):
    """
    return extension of results from plink
    This is handy when the extension deviates from ".assoc"
    for instance when --covar option is used
    """
       
    shortened_list = []
    for plink_out in list_of_files:
        if not plink_out.endswith("adjusted"):
            if not plink_out.endswith("log"):
                shortened_list.append(plink_out)
    
    regex = re.compile("P[0-9]+\.(?P<extension>[A-z\.]+)$", re.MULTILINE)
    extensions = set(regex.findall("\n".join(shortened_list)))
    
    if len(extensions) == 1:
        extension_str = extensions.pop()
        
    elif len(extensions) == 0:
        regex = re.compile("\.(?P<extension>[A-z\.]+)$", re.MULTILINE)
        extensions = set(regex.findall("\n".join(list_of_files)))
        extension_str = extensions.pop()
        
    else:
        log.info(str(extensions))
        raise IOError(len(extensions), str(extensions) + str(list_of_files))

    return(extension_str)


def check_bim_bed_fam(path):
    """
    check fam , bim and bed file is present
    """

    if path[-4:] in [".fam", ".bim", ".bed"]:     #check path ends with ".fam",".bim",".bed" and remove this extension
        path = path[:-4]

    getfile_handle(path + ".fam" , "This is --bfile variable", False)
    getfile_handle(path + ".bed" , "This is --bfile variable", False)
    getfile_handle(path + ".bim" , "This is --bfile variable", False)

    return path


def getfile_handle(path, message="", verbose=True):
    """
    return a filehandle of path or quit
    """
    
    if not path:
        log.info(message + '\nPath not defined for ')
        if verbose:
            raise IOError()
        else:
            endtime = time.asctime(time.localtime(time.time()))
            sys.exit("\nTerminated analysis: " + endtime + "\n")
    try:
        filehandle = open(path, 'rU')
        if(os.path.getsize(path) == 0):
            log.info("\nWarning: file size of " + path + " is zero")
        return filehandle
    
    except IOError, error:
        log.info("\n" +  "\nCannot find " + path) 
        if verbose:
            raise error
        else:
            endtime = time.asctime(time.localtime(time.time()))
            sys.exit("\nTerminated analysis: " + endtime + "\n")


def remove_plink_output(resultfile):
    """
    removes all files created by plink. 
    should function the same as bash "rm resultsfile*":
    for instance "rm /tmp/rndvar*"
    """
    files_create_by_plink = glob.glob(resultfile + "*")
    for file_to_delete in files_create_by_plink:
        os.remove(file_to_delete)


def average(numbers):
    """
    calculate average (mean) of a list of numbers
    """
    numbers = array.array('d', numbers)

    return (sum(numbers) / len(numbers))


def variance(numbers):
    """
     calculate sample variance of a list of numbers.
     equal to R-projects var() implementation
    """
    if not isinstance(numbers, list):
        raise AssertionError

    numbers = array.array('d', numbers)

    if(len(numbers) < 2):
        return "NA"
    else:
        numbers = [float(x) for x in numbers]
        mean = average(numbers)
        diffrence_from_mean = [a - b for a, b in zip(numbers, [mean] * len(numbers))]
        sum_of_diff = sum([a ** b for a, b in zip(diffrence_from_mean, [2] * len(numbers))])
        var = ((1 / (float(len(numbers)) - 1)) * sum_of_diff)

    return var


def write_text_to_tempfile(text_file):
    """
    writes a file to temporary file (in linux to /tmp/) and returns the filename
    """
    filename = tempfile.mkstemp()
    file_handle = open(filename[1], "w")
    # Write all the lines at once:
    file_handle.writelines(text_file)
    file_handle.close()
    return filename[1]


def extract_pheno_from_pheno_file(pheno_path, pheno_available):
    """
    extract pheno type from a pheno_file and include only the columns which are present in pheno_available
    as pheno_perm_clean and the start of the lines as array 
     
    """
    pheno_arrays = read_pheno_file(pheno_path)
    pheno_file_list = pheno_arrays[0] #this are the IDs
    pheno_available_as_index = [x - 1 for x in pheno_available]
    pheno_available_as_index.sort()
    pheno_perm = [x.split("\t") for x in pheno_arrays[1]]

    if(len(pheno_perm[0]) != len(pheno_available_as_index)):
        log.info("WARNING: Amount of sumlog files (" + str(len(pheno_available_as_index)) \
              + ") is not the same as amount of phenotypes (" + \
              str(len(pheno_perm[0])) + \
              ") defined in phenotype file")

    pheno_perm_list_clean = [[y[w] for w in pheno_available_as_index] for y in pheno_perm]
    pheno_perm_clean = ["\t".join(q) for q in pheno_perm_list_clean]

    return  pheno_perm_clean, pheno_file_list


def extract_pheno_from_fam_file(path_fam_bed_bim):
    """
    extract pheno type from a fam (sixth column and further ) as pheno_perm_clean and 
    the start of the lines(first two columns) as array 
    """
    fam_fh = getfile_handle(path_fam_bed_bim + ".fam")
 
    pheno_perm_clean, pheno_file_list = split_fh_in_columns(fam_fh, 2, 5)
    return pheno_perm_clean, pheno_file_list


def split_fh_in_columns(covar_fh, split, split_2=None):
    """
    function to split a filehandle in a column orientation manner
    returns a tuple of 2 lists of strings
    The character to split on is tab
    the split parameter is a integer to pinpoint the split
    columns in the middle can be omited by using  spit_2     
    """
    if split_2 == None:
        split_2 = split

    covar_tuple = [(t[:split], t[split_2:]) for t in [x.strip().split() for x in covar_fh]]
    covar_perm_clean = [" ".join(x[1]) for x in covar_tuple]
    covar_file_list = [" ".join(x[0]) for x in covar_tuple]

    return covar_perm_clean, covar_file_list


def extract_pheno_from_covar_file(covar_file):
    """
    extract pheno type from a covar as pheno_perm_clean and 
    the start of the lines(first two columns) as array 
    """
    covar_fh = getfile_handle(covar_file)
    covar_perm_clean, covar_file_list = split_fh_in_columns(covar_fh, 2)
    return covar_perm_clean, covar_file_list



def make_shuffled_phenotype(seeds, pheno_perm_clean, start_of_lines):
    """
    """
    for i in range(len(seeds)):
        random.seed(seeds[i])
        random.shuffle(pheno_perm_clean)
        for j in range(len(start_of_lines)):
            start_of_lines[j] = start_of_lines[j] + " " + str(pheno_perm_clean[j])

    perm_text_file = ""
    
    for i in range(len(start_of_lines)):
        perm_text_file = perm_text_file + start_of_lines[i] + "\n"

    filename = write_text_to_tempfile(perm_text_file)
    return filename


def create_pheno_permutation_file(plink, seeds):
    """
    create permuted pheno file where length of the list seeds is the amount of permutations.
    pheno available is a list of colum numbers that must be permuted (same as in clusterresults)
    if no pheno file is present the 6th column of the map file will be used
    """
    
    pheno_path = plink.pheno_file
    path_fam_bed_bim = plink.bfile
    pheno_available = plink.phenotypes_present
   
    if (not pheno_path == None):
        pheno_perm_clean, start_of_lines = extract_pheno_from_pheno_file(pheno_path, pheno_available)         #extract from pheno file
    else:
        # get 6th column(and latter) as clean phenotype
        pheno_perm_clean, start_of_lines = extract_pheno_from_fam_file(path_fam_bed_bim)
    
    filename = make_shuffled_phenotype(seeds, pheno_perm_clean, start_of_lines)

    return(filename)


def create_covar_permutation_file(plink, seeds):
    """
    """
    covar_perm_clean, start_of_lines = extract_pheno_from_covar_file(plink.covar_file)
    filename = make_shuffled_phenotype(seeds, covar_perm_clean, start_of_lines)
    return (filename)
    

def save_text_file(filename, text):
    """ saves a string to a filename
    If the file already exists the text will be printed and is not saved
    """
    file_handle = open(filename, "w")
    file_handle.write(text)
    file_handle.close()


def read_pheno_file(path):
    """
    read a pheno file from path and return first 2 columns as list
    and the remaining columns as list
    """

    file_handle = getfile_handle(path)

    start_of_lines = []  #read the cluster of genes file
    phenotypes = []
    
    for line in file_handle:
        splitted = line.split()
        start_of_lines.append(splitted[0] + " " + splitted[1])
        phenotypes.append("\t".join([str(q) for q in splitted[2: ]]))

    return([start_of_lines, phenotypes])


def map_snp_to_geneset(geneset, gene_group=False):
    """
    read the mapping from snp id to gene and clusternumber
    the output is a dict with primary key the snp name
    each value in the primary dict consist out of a dict with a gene (key="g")
     and a cluster number (key="c") 
    """
    
    if len(geneset) == 0:
        sys.exit("no value for geneset file selected")
    
    file_handle = getfile_handle(geneset)
    snptogeneset = {}
    #read the cluster of genes file
    for line in file_handle:
        splitted = line.split()
        #assign variables to descriptive variable name
        if (len(splitted) > 2):
            snp = splitted[0]
            gene = splitted[1]
            cluster = str(splitted[2])
        elif(len(splitted) == 2):
            snp = splitted[0]
            gene = ''
            cluster = splitted[1]
        else:
            if len(line.strip()) == 0: # if blank lines
                pass
            else:
                log.info("geneset file has only one 1 column: this should be at least three with on :\
                \ncolumn 1 the SNP \ncolumn 2 the geneid\
                \ncolumn 3 geneset (without spaces in the name)")
                log.info("the following line is wrong: " + line)
                sys.exit(1)

        if (not snp in snptogeneset):
            snptogeneset[snp] = []
        if not gene_group:
            snptogeneset[snp].append({"g":gene, "c":cluster})
        else:
            snptogeneset[snp].append({"g":gene, "c":gene})

    return  snptogeneset


def map_geneset_to_snp(geneset):
    """
    read the mapping from snp id to gene and geneset number
    the output is a dict with primary key the snp name
    each value in the primary dict consist out of a dict with a gene (key="g")
    and a geneset number (key="c") 
    """
    if len(geneset) == 0:
        sys.exit("no value for geneset file selected")
    file_handle = getfile_handle(geneset)
    snptogeneset = {}
    
    for line in file_handle:
        splited = line.split()
        #assign variables to descriptive variable name
        snp = splited[0]
        gene = splited[1]
        cluster = splited[2]

        if (not cluster in snptogeneset):
            snptogeneset[cluster] = []
        snptogeneset[cluster].append({"g":gene, "s":snp})

    return  snptogeneset


def negsumlog10(numbers):
    """
    calculate negsumlog10 on a list of numbers
    """
    numbers = array.array('d', numbers)
    neg_sum_log10_value = sum([ math.log10(x) for x in numbers])
    return abs(neg_sum_log10_value)

def alpha_sort(key):
    parts = re.split('(\d*\.\d+|\d+)', key)
    return tuple((e.swapcase() if i % 2 == 0 else float(e))
            for i, e in enumerate(parts))


def get_starttime():
    """
    get system start time
    """
    starttime = time.asctime(time.localtime(time.time()))
    return("Analysis started: " + starttime)



def get_endtime():
    """
    get system end time
    """
    endtime = time.asctime(time.localtime(time.time()))
    return("\nFinished analysis: " + endtime + "\n")



def get_terminated_time():
    """
    get system end time
    """
    termtime = time.asctime(time.localtime(time.time()))
    return("\nAnalysis terminated: " + termtime + "\n")

