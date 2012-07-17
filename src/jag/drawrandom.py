"""
Make random draws

@author: maartenk, esther lips
"""

import jag.common as common
import sys
import random
import logging as log
import os.path
import time


def _get_terminated_time():
    """
    get system end time
    """
    termtime = time.asctime(time.localtime(time.time()))
    return("\nAnalysis terminated: " + termtime + "\n")


class DrawRandom(object):
    """
    Class to draw random genes or snps
    """

    def __init__(self, genemapping, snpmapping, inoutput):
        self.genemapping = genemapping
        self.snpmapping = snpmapping
        self.inoutput = inoutput
        self.accessed = None
        self.exclude = True


    def setseed(self, seed):
        """
        setting seed to a value: this function is used for testing and debugging
        """
        random.seed(seed)
        log.info("setting seed to value: " + str(seed))


class Genes(object):
    """
    class for drawing SNP's out of genes
    """
    def __init__(self, drawrandom):
        self.drawrandom = drawrandom
        self.exclude_genes = {}


    def _gene2snpmapping(self):
        """
        create a mapping from a gene to SNP mapping file
        The key a is a gene and the value is a list of SNP's
        """
        
        allsnp_and_genes_fh = common.getfile_handle(self.drawrandom.genemapping)
        gene2snp_mapping = {}
        for text in allsnp_and_genes_fh:
            text_array = text.strip().split("\t")
            
            if len(text.strip()) != 0:  # if line is not empty
                try:
                    gene2snp_mapping[text_array[1]].append(text_array[0])
                except KeyError:
                    gene2snp_mapping[text_array[1]] = [text_array[0]]

        return gene2snp_mapping


    def genes(self, geneset, amount):
        """
        Draw random genes, where set is the geneset name selected
         and amount the number of new sets to be formed
        """
        geneset_to_snp_mapping = common.map_geneset_to_snp(self.drawrandom.snpmapping)
        if (geneset not in geneset_to_snp_mapping):
            log.info("\n" + geneset + " is not known as a geneset name. Please check your data for correct geneset name.")
            log.info(_get_terminated_time())
            sys.exit()
        #load number of  genes in geneset by the mapping file
        self.exclude_genes = set([ x["g"] for x in geneset_to_snp_mapping[geneset]])

        ngenes = len(self.exclude_genes)
        gene2snp = self._gene2snpmapping()  #create gene to all SNP's mapping
        all_genes = set(gene2snp.keys())
        
        if(len(all_genes) <= (ngenes*amount)):
            max = len(all_genes)/ngenes
            log.info("\nYour sample of " + str(len(all_genes)) + " genes is too small to draw " + str(amount) + " random genesets of size " + str(ngenes))
            log.info("\nMaximum number of independent gene-sets to draw (" + str(len(all_genes))+"/"+ str(ngenes) + "): " + str(max))
            log.info(_get_terminated_time())
            sys.exit()
        
        
        log.info("\nSize of gene pool (--pool): " + str(len(all_genes)) + " unique genes")
        log.info("\nThe geneset " + geneset + " contains " + str(ngenes) + " genes" )
        log.info("\nMinimum size of gene pool needed to draw " + str(amount) + " gene-sets (" + str(ngenes) + "*" +  str(amount) + "): " + str(ngenes*amount))
   
        #exclude genes that are in the selected geneset
        if (self.drawrandom.exclude):
            all_genes = all_genes.difference(set(self.exclude_genes))

        output_text = ""
        
        for new_set_number in xrange(1, amount + 1):

            setname = "Draw_" + str(new_set_number)
            sampled_keys = random.sample(all_genes, ngenes)

            for gene_as_key in sampled_keys:
                output_text += '\n'.join([str(snp) + "\t" + str(gene_as_key) + \
                        "\t" + setname for snp in gene2snp[gene_as_key]])
                output_text += '\n'

        incl_or_excl = "excl" if self.drawrandom.exclude else "incl"      #create correct filename
        filename = "draws_ngenes.set.annot"
        out = self.drawrandom.inoutput.save_text_to_filename(filename, output_text)     #save file
        log.info("\nSaved random draws on number of genes as " + out)
        
        return(output_text)


class Snp(object):
    """
    class for drawing SNP's
    """
    def __init__(self, drawrandom):
        self.drawrandom = drawrandom
        self.inregion = "all"
        self.exclude_snps = {}


    def _create_indep_snp_file(self, plink):
        """
        Create a list of independent snps.
        plink: a plink object.(only needs a bfile)
        
        """
        log.info("\nPerforming LD based pruning...")
        plink.set_plink_arguments("--bfile " + plink.bfile + " --indep-pairwise 200 5 0.25")
        resultfile = plink.run_plink()
        prune_file = resultfile + ".prune.in"
     
        snp2gene_mapping = {}   #create genemapping mapping
        file_handle = common.getfile_handle(self.drawrandom.genemapping)
        
        for text in file_handle:
            text_array = text.strip().split()
            if len(text.strip()) != 0:
                try:
                    snp2gene_mapping[text_array[0]] += "," + text_array[1]
                except KeyError:
                    snp2gene_mapping[text_array[0]] = text_array[1]
                        
        # use snp2gene_mapping to map the snp's to genes
        outfile_text = ""
        file_handle = common.getfile_handle(prune_file)

        for text in file_handle:
            rs_number = text.strip()
            try:
                outfile_text += rs_number + '\t' + str(snp2gene_mapping[rs_number]) + "\n"  #add mapping
            except KeyError:
                outfile_text += str(rs_number) + "\t" + " " + "\n"  
        
        outfile = self.drawrandom.inoutput.save_text_to_filename("prune.in", outfile_text)  #save mapped prune file to prunefile
        log.info("Saved pruned SNP file as " + outfile)
        
        return (outfile)


    def _get_random_snp(self, allsnp_and_genes, amount_allsnp_and_genes):
        """
        draw SNPs or genes without replacement (drawn items in self.accessed)
        Also checks optional if gene or snp is in original list group mapping and exclude
        this gene/snp to be drawn This is regulated by self.exclude (True of False)
        
        allsnp_and_genes: list of raw readlines in form "RS12345\t12345\n"
        amount_allsnp_and_genes: length of amount_allsnp_and_genes
        """
        
        rand_acces = random.randint(0, amount_allsnp_and_genes)
        
        if (rand_acces not in self.drawrandom.accessed):
            self.drawrandom.accessed.add(rand_acces)
            rsnumber, genename = [a.strip() for a in allsnp_and_genes[rand_acces].split('\t')]
            
            if(self.inregion == "all"):
                pass
            
            elif(self.inregion == "genic"):
                if(genename == ""):
                    rsnumber, genename = self._get_random_snp(allsnp_and_genes, amount_allsnp_and_genes)
                    
            elif(self.inregion == "intergenic"):
                if(genename != ""):
                    rsnumber, genename = self._get_random_snp(allsnp_and_genes, amount_allsnp_and_genes)

            if(self.drawrandom.exclude):    #check gene is in the geneset in snp mapping otherwise draw new
                if(rsnumber in self.exclude_snps): #allsnp_and_genes.pop(rand_acces)
                    rsnumber, genename = self._get_random_snp(allsnp_and_genes, amount_allsnp_and_genes)

        else:
            rsnumber, genename = self._get_random_snp(allsnp_and_genes, amount_allsnp_and_genes)

        return rsnumber, genename


    def snp(self, geneset, amount, empp_file, plink, out):
        """
        draw random snp, on number of neff snps in the .empp file, from an independent snp file (.prune.in). If this file is not 
        present, this file will be created within _create_indep_snp_file function.
        
        """
        
        geneset_to_snp_mapping = common.map_geneset_to_snp(self.drawrandom.snpmapping)
        
        if (geneset not in geneset_to_snp_mapping):
            log.info("\n" + geneset + " is not known as a geneset name. Please check your data for correct geneset name.")
            log.info(_get_terminated_time())

            sys.exit()

        self.exclude_snps = set([ x["s"] for x in geneset_to_snp_mapping[geneset]]) #drawsnps
        
        try:
            n_snp = _read_empp_results(empp_file)[geneset]

        except ValueError:
            log.info("geneset to select is not found in empp file")
            log.info(_get_terminated_time())
            
            sys.exit()
        
        prunedin_file = str(os.path.abspath(os.path.curdir)) + "/" + out + ".prune.in" #get pruned.in file
        prunedin = os.path.exists(prunedin_file)
        
        if (prunedin == False): #create prune.in from file, if it does not excist
            prunedin = self._create_indep_snp_file(plink)
            
        else:
            log.info ("\nUsing the pruned SNP set from " + prunedin_file)
            prunedin = prunedin_file

        allsnp_and_genes = common.getfile_handle(prunedin).readlines() #pylint: disable=E1103
        
        length_genic = 0
        length_nongenic = 0
        
        for line in allsnp_and_genes:
            line = line.rsplit('\t')
            gene = line[1].rsplit()
                     
            if len(gene) is not 0:
                length_genic += 1
            else:
                length_nongenic += 1
                                     
        amount_allsnp_and_genes = len(allsnp_and_genes) - 1
              
        #random_snp_text = "RS#\tGeneID\tDraw_#\n"
        random_snp_text = ""
        
        for new_group_number in xrange(1, amount + 1):
            self.drawrandom.accessed = set()
            setname = "Draw_" + str(new_group_number)

            count = 0
            while(count < n_snp):
                random_snp_text += ("\t".join(self._get_random_snp(allsnp_and_genes, amount_allsnp_and_genes)))
                random_snp_text += ("\t" + setname + "\n")
                count = count + 1
                
        inregion_text = "unknown_in_region"
        
        if(self.inregion == "genic"):
            inregion_text = "genic"
            
            if amount*n_snp > length_genic:
                log.info("\nWarning: pool of " + str(length_genic) + " SNPs located within genes is to small to draw " + \
                 str(amount) + " x " + str(n_snp) + " independent nEff SNPs!")
                log.info(_get_terminated_time())
                sys.exit()
            else:               
                log.info("\nDrawing " + str(amount) + " x " + str(n_snp) + " nEff SNPs from a pool of " + \
                 str(length_genic) + " SNPs located within genes")
            
        elif(self.inregion == "intergenic"):
            inregion_text = "intergenic"
            
            if amount*n_snp > length_nongenic:
                log.info("\nWarning: pool of " + str(length_nongenic) + " SNPs located outside genes is to small to draw " + \
                 str(amount) + " x " + str(n_snp) + " independent nEff SNPs!")
                log.info(_get_terminated_time())
                sys.exit()
            else:
                log.info("\nDrawing " + str(amount) + " x " + str(n_snp) + " nEff SNPs from a pool of " + \
                 str(length_nongenic) + " SNPs located outside genes")
            
        elif(self.inregion == "all"):
            inregion_text = "all"
            
            if amount*n_snp > len(allsnp_and_genes):
                log.info("\nWarning: pool of " + str(len(allsnp_and_genes)) + " SNPs located in- and outside genes is to small to draw " + \
                 str(amount) + " x " + str(n_snp) + " independent nEff SNPs!")
                log.info(_get_terminated_time())
                sys.exit()
            else:
                log.info("\nDrawing " + str(amount) + " x " + str(n_snp) + " nEff SNPs from a pool of " + \
                 str(len(allsnp_and_genes)) + " SNPs located in- and outside genes")
            
        incl_or_excl = "unknown"
        
        if (self.drawrandom.exclude is False):
            incl_or_excl = "incl"
            
        elif (self.drawrandom.exclude is True):
            incl_or_excl = "excl"
        
        filename = "draws_neff_" + inregion_text + ".set.annot"
        out = self.drawrandom.inoutput.save_text_to_filename(filename, random_snp_text)     #save random snps file
        log.info("\nSaved random draws on number of effective number of SNPS as " + out)
        return(filename)


def _read_empp_results(file_name):
    """
    Read a empp file and get nEff back as a dict with the genesetname as key
    """

    file_handle = common.getfile_handle(file_name)
    header_raw = file_handle.readline()
    header = header_raw.strip().split("\t")

    try:
        neff_column = header.index("nEff")
    except ValueError:
        sys.exit("nEff column not found in header of " + file_name)

    geneset_neff_mapping = {}
    
    for line in file_handle:
        splitted_line = line.strip().split("\t")
        geneset_neff_mapping[str(splitted_line[0])] = int(splitted_line[neff_column])
    
    return (geneset_neff_mapping)


