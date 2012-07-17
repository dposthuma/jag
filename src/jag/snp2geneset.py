'''
Created on May 9, 2011

@author: maartenk, esther lips
'''

import jag.common as common
import bisect
from operator import itemgetter
import logging as log


def _read_snp_file(snp_data_file, all_chromosomes):
    """
    Read snp_data file and skip the chromosomes that are not
    present in the list all_chromosomes
    """
    chrsnpmapping = {}
    for chr_name in all_chromosomes:
        chrsnpmapping[chr_name] = []

    snp_file_handle = common.getfile_handle(snp_data_file)
    old_chr = ""

    for line in snp_file_handle:
        try:
            snp_data = line.strip().split("\t")
            current_chr = snp_data[1]
 
            if not current_chr in all_chromosomes:                
                if not old_chr == current_chr:
                    log.info("There are no genes on chromosome " + current_chr)
              
                old_chr = current_chr
                
            else:
                if not old_chr == current_chr:
                    log.info("Mapping SNPs to genes on chromosome " + str(current_chr))
                    old_chr = current_chr

                #save snp
                rsnumber = snp_data[0]
                try:
                    chr_location = int(snp_data[2])
                    chrsnpmapping[current_chr].append((chr_location, rsnumber))

                except ValueError:
                    log.info("Ignoring following line in SNP file: " + "\t".join(snp_data))
                    
        except IndexError:
            log.info("indexerror" + str(line))
            
        except ValueError:
            log.info (line)
            
    log.info("")
    
    return(chrsnpmapping)


def find_snps(chrsnpmapping, current_chromosome, sorted_keys, gene, leftside, rightside):
    """
    Map snps to genes
    """
    #searching for snps
    left_index = bisect.bisect_left(sorted_keys, leftside)
    right_index = bisect.bisect_right(sorted_keys, rightside)

    if (left_index == 0 and right_index == 0):
        pass

    snps_raw = chrsnpmapping[current_chromosome][left_index:right_index]
  
    return snps_raw


def _lookupgenes(start_ext_base, end_ext_base, chrsnpmapping, all_genes_from_all_groups):
    """
    lookup all genes in the list all_genes_from_all_groups.
    The SNPs must be saved in chrsnpmapping  which is made by _read_snp_file()
    The upstream and downstream SNPs of a gene which is included in the mapping
    is specified by  start_ext_kb and end_ext_kb and should be in base
    """
    current_chromosome = None
    sorted_keys = []
    
    for gene in all_genes_from_all_groups:
        if not current_chromosome == gene["chr_number"]:
            current_chromosome = gene["chr_number"] #order chromosome
            chrsnpmapping[current_chromosome].sort(key=lambda r:r[0]) #create search on key list
            sorted_keys = [r[0] for r in chrsnpmapping[current_chromosome]]
        #calculate start and end of gene
        leftside = 0
        rightside = 0
        
        if gene["chr_strand"] == "+":
            leftside = gene["chrstart"] - start_ext_base
            rightside = gene["chrend"] + end_ext_base
        
        elif gene["chr_strand"] == "-":
            rightside = gene["chrend"] + start_ext_base
            leftside = gene["chrstart"] - end_ext_base
        
        else:
            log.info("Invalid strand position" + str(gene["chr_strand"]))
        
        snps_raw = find_snps(chrsnpmapping, current_chromosome, sorted_keys, gene, leftside, rightside)
        gene["snp_list"] = set([rs_tuple[1] for rs_tuple in snps_raw])
        
    return all_genes_from_all_groups


def _read_group_gene_mapping(infile):
    """
    read file with genes for gene mapping
    """
    gene_group_mapping = [] #read gene file
    file_handle = common.getfile_handle(infile)
    
    for line in file_handle: 
        gene_group_mapping.extend([line.strip().split("\t")])
        
    return gene_group_mapping


def unique(seq): 
    """
    get uniques
    """ 
    # order preserving 
    checked = [] 
    for e in seq: 
        if e not in checked: 
            checked.append(e) 
    return checked


def _create_gene2group_mapping(gene_group_mapping):
    """
    map genes per gene-set
    """
    group2gene = {}
    # get unique group members
    un_gene_group_mapping = unique(gene_group_mapping)
    
    for line in un_gene_group_mapping:
        try:
            group2gene[line[1]].append(line[0])
        except KeyError:
            group2gene[line[1]] = [line[0]]

    return group2gene


def _read_gene_boundaries(gene_boundaries_file, genes_in_groups):
    """
    # read gene boundaries
    """
    gb_file = common.getfile_handle(gene_boundaries_file)
    genes_with_boundaries = []
    #header = gb_file.readline() 

    for line in gb_file:
        gene_data = line.strip().split("\t")
     
        if gene_data[0] in genes_in_groups:
            genes_in_groups.remove(gene_data[0]) #save gene information
            chr_number = gene_data[1]
            chrstart = int(gene_data[2])
            chrend = int(gene_data[3])
            chr_strand = gene_data[4]
            symbol = gene_data[1]
            gene_info = {"geneid":gene_data[0],
                        "chr_number":chr_number,
                        "chrstart":chrstart,
                        "chrend":chrend,
                        "chr_strand":chr_strand,
                        "symbol":symbol,
                        "snp_list":[]
                        }
            genes_with_boundaries.append(gene_info)
    
    genes_without_boundaries = genes_in_groups
    return(genes_with_boundaries, genes_without_boundaries)


def _remove_genes_without_boundaries_from_groups(group2gene, non_found_genes):
    """
    """
    nbSet = [] # set of genes without gene boundaries
    
    if not len(non_found_genes) == 0:
        log.info("No gene boundaries are found for:")
        for gene in non_found_genes:
            for genes in group2gene.items():
                if gene in genes[1]:
                    nbSet.append(gene)
                    genes[1].remove(gene)
        
    unique = set(nbSet)
    for item in unique:
        log.info("geneID " + item)
    
    log.info("")
    return group2gene


def _map_lookedup_snps_to_group(all_genes_from_all_groups, group2gene):
    """
    """
    geneid2gene = {}
    snp2all_gene = {}
    genelist = []
    
    for gene in all_genes_from_all_groups:
        geneid2gene[gene["geneid"]] = gene
        info = [gene["geneid"], gene["snp_list"]]
        genelist.append(info)
    
    snp2all_gene['CODING_GENES'] = genelist 
    
    group2snp = {}
    
    for (group, genes) in group2gene.iteritems():
        try:
            allsnps_per_group = [(geneid, geneid2gene[geneid]["snp_list"]) for geneid in genes]
            group2snp[group] = allsnps_per_group

        except KeyError:
            log.info("")
            log.info("Cannot not find SNPs for geneID " + str(genes))
        
    return snp2all_gene, group2snp


def format_mapping(group2snp):
    output_text = ""

    for group, snps_tupple in group2snp.iteritems():
        for gene_and_snps in snps_tupple:
            snps = gene_and_snps[1]
            gene = gene_and_snps[0]

            if(len(snps) is not 0): # to avoid white lines in output
                output_text += '\n'.join([snp + "\t" + gene + "\t" + group for snp in snps])
                output_text += '\n'
    return output_text


def snp2geneset(infile, start_ext_kb, end_ext_kb, gene_boundaries_file, snp_data_file):
    """
    MAIN annotate snps to genes
    """
    # calculate start extra in base pair instead in kilo base pair
    start_ext_kb = float(start_ext_kb) * 1000
    end_ext_kb = float(end_ext_kb) * 1000

    # read gene-set file
    gene_group_mapping = _read_group_gene_mapping(infile)    
    group2gene = _create_gene2group_mapping(gene_group_mapping) #create a group to gene mapping
    genes_in_groups = set([group[0] for group in gene_group_mapping]) #find all unique genes
       
    # read data for all genes
    all_genes_to_map = _read_group_gene_mapping(gene_boundaries_file) # read in datafile for all genes    
    all_genes = set([gene[0] for gene in all_genes_to_map])
       
    log.info("Found " + str(len(all_genes)) + " unique genes in " + gene_boundaries_file ) 
    log.info("Found " + str(len(genes_in_groups)) + " unique genes in " + infile + "\n") 
    
    # remove genes from gene-sets that do not have gene boundary info
    genes_with_boundaries, genes_without_boundaries = _read_gene_boundaries(gene_boundaries_file, genes_in_groups)
    group2gene = _remove_genes_without_boundaries_from_groups(group2gene, genes_without_boundaries)
     
    # remove genes from all genes that do not have gene boundary info
    all_genes_with_boundaries, all_genes_without_boundaries = _read_gene_boundaries(gene_boundaries_file, all_genes)
    all_genes_with_boundaries = sorted(all_genes_with_boundaries, key=itemgetter('chr_number'))
    all_chromosomes = set([gene["chr_number"] for gene in all_genes_with_boundaries])

    chrsnpmapping = _read_snp_file(snp_data_file, all_chromosomes)  # read snp file
    log.info("Saving results...")
     
    all_genes_with_snps = _lookupgenes(start_ext_kb, end_ext_kb, chrsnpmapping, all_genes_with_boundaries)

    snp2all_genes,group2snp = _map_lookedup_snps_to_group(all_genes_with_snps, group2gene)
    
    return(snp2all_genes, group2snp)


