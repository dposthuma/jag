'''
Created on Dec 14, 2010

@author: maartenk, esther lips

Stores, reads and write results and calculate concatenated test statistics
'''

import jag.common as common
import sys
import logging as log



class Clusterresults(object):
    '''
    Stores and reads and write results and calculate concatenated scores
    '''

    def __init__(self):
        self.cluster_results = {}
        self.column_number = None
        self.total_columns = None
        self.pheno_path = None
        self.geneset_path = None
        self.bfile_path = None

        self.clean_results = {}


    def add_snp_to_cluster_all_info(self, p_value, snp_id, gene_and_cluster_list):
        """
        add a P value found in PLINK output to the cluster in gene_and_cluster_list
        """
        for cluster in gene_and_cluster_list:
            if  not cluster["c"] in self.cluster_results:
                self.cluster_results[cluster["c"]] = {"genes":[], "snp":[], "p-value":[]}
                
            self.cluster_results[cluster["c"]]["genes"].append(cluster["g"])
            self.cluster_results[cluster["c"]]["snp"].append(snp_id)
            
            if(p_value == "NA"):
                p_value = 1
               
            self.cluster_results[cluster["c"]]["p-value"].append(float(p_value))
                
            if (p_value == "P"):
                log.info("The set (--set) file should not contain a header")
                sys.exit("\nTerminated analysis: " + endtime + "\n")
 


    def write_all_pvalues(self):
        """
        Write all P values and return file containing P-values each on a new line.
        """
        
        message = []
        for group  in self.cluster_results.iterkeys():
            pvalues = self.cluster_results[group]["p-value"]
            message.extend([group + "\t" + str(pvalue) for pvalue in pvalues])

        filename = common.write_text_to_tempfile("\n".join(message))

        return(filename)


    def map_p_values_from_assoc_file(self, snptogroup, assoc_file, adjusted=False):
        """
        extracts P-values from a assoc file and adds them based on the 
        snptogroup mapping to the right gene-set
        
        """
        file_handle = common.getfile_handle(assoc_file)
        
        for line in file_handle:
            splitted = line.split()
            snp_id = splitted[1]
            
            if adjusted:
                p_value = splitted[3]
                                
            else:
                p_value = splitted[8]
                                
            if (snp_id in snptogroup):
                self.add_snp_to_cluster_all_info(p_value, snp_id, snptogroup[snp_id])


    def createcleanresults(self):
        """
        summarize cluster_results and save it in self.clean_results
        """
       

        for cluster in self.cluster_results.iterkeys():
            current_cluster = self.cluster_results[cluster]
            n_snps = len(set(current_cluster["snp"]))
            n_genes = len(set((current_cluster["genes"])))
            nsumlog10 = common.negsumlog10(current_cluster["p-value"])

            self.clean_results[cluster] = {"n_snps":n_snps,
                                        "n_genes":n_genes,
                                        "sumneglog10":nsumlog10}
        
        return(self.clean_results)


    def get_sum_neglog_10_score(self, group):
        """
        retrieve sum neglog 10 value if not found return -1
        """
        score = -1
        if group in self.clean_results:
            score = self.clean_results[group]["sumneglog10"]
        else:
            sys.exit("could not find -log10 for set with name " + str(group))
        return score


    def get_n_snps(self, group):
        """
        retrieve amount of snps if not found return -1

        """
        n_snp = -1
        if group in self.clean_results:
            n_snp = self.clean_results[group]["n_snps"]
        else:
            sys.exit("could not find score for set with name " + str(group))
        return n_snp


    def get_n_genes_of_group(self, group):
        """
        get number of genes that is in the group
        """
        n_genes = -1
        if group in self.clean_results:
            n_genes = self.clean_results[group]["n_genes"]
        else:
            sys.exit("could not find score for set with name " + str(group))
        return n_genes


    def format_results(self):
        """
        format the results in a way that a user can read AND read_formated_results can read
        """
        
        formated_text = "phenotype #\t" + str(self.column_number)
        formated_text += "\n# of phenotypes in pheno file\t" + str(self.total_columns)
        formated_text += "\npheno file path\t" + str(self.pheno_path)
        formated_text += "\ngroup file path\t" + str(self.geneset_path)
        formated_text += "\nbfile path\t" + str(self.bfile_path)
        formated_text += "\nGroup\tnSNPs\tnGENES\tsum(-log10(p))\n"
        
        keys = self.clean_results.keys()
        keys.sort(key=common.alpha_sort) 
        
        for cluster in keys:
            current_cluster = self.clean_results[cluster]
            formated_text += str(cluster)
            formated_text += "\t" + str(current_cluster["n_snps"])
            formated_text += "\t" + str(current_cluster["n_genes"])
            formated_text += "\t" + str(current_cluster["sumneglog10"]) + "\n"
            
        return (formated_text)


    def read_formated_results(self, file_handle):
        """
        import the results saved by format_results() and populate this class
        """

        text_as_list = file_handle.readlines()

        if len(text_as_list) == 0:
            raise AssertionError
        self.column_number = text_as_list[0].strip().split("\t")[1]
        self.total_columns = text_as_list[1].strip().split("\t")[1]
        self.pheno_path = text_as_list[2].strip().split("\t")[1]
        self.geneset_path = text_as_list[3].strip().split("\t")[1]
        self.bfile_path = text_as_list[4].strip().split("\t")[1]

        text_clean = [c.split("\t") for c in[ l.strip() for l in text_as_list[5:]]]
        
        for i in range(1, len(text_clean)):
            self.clean_results[text_clean[i][0]] = {"n_snps":text_clean[i][1],
                                        "n_genes":text_clean[i][2],
                                        "sumneglog10":text_clean[i][3]}



