#!/usr/bin/env python

'''
Created on Dec 13, 2010

@author: maartenk, estherlips

main class where program is called with from command line
'''

import getopt
import sys
import os
import logging as log
import time
import re

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/src/")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../") #for development reasons

import jag.common as common

def _get_starttime():
    """
    get system start time
    """
    starttime = time.asctime(time.localtime(time.time()))
    return("Analysis started: " + starttime)


def _get_endtime():
    """
    get system end time
    """
    endtime = time.asctime(time.localtime(time.time()))
    return("\nFinished analysis: " + endtime + "\n")
    
    
def _get_terminated_time():
    """
    get system end time
    """
    termtime = time.asctime(time.localtime(time.time()))
    return("\nAnalysis terminated: " + termtime + "\n")
    
    
def _get_header_of_program():
    """
    print the header of the program
    """
    
    return("""
    ##############################################################
    #                                                            #
    #     JAG  - Joint Association of Genetic variants - V1.0    #   
    #            2012,  Esther Lips & Maarten Kooyman            #
    #              GNU General Public License, v2                # 
    #                                                            #   
    #                Complex Trait Genetics Lab                  #
    #               http://ctglab.nl/software/jag/               #
    #                                                            #
    ##############################################################
    """)


def _usage():
    """
    print help for using this program

    """
    return("""
    === SNP2GENE ANNOTATION ===
    --snp2gene         file with genes to map
    --up               size upstream region in kb
    --down             size downstream region in kb
    --gene_loc        file with chromosomal locations of genes 
    --snp_loc         file with chromosomal location of snps
     
    === SELF-CONTAINED TEST ===
    --bfile    -b      path of fam, bim and bed file without file extension
    --set      -g      file with snp annotated genesets 
    --perm     -m      number of permutations
    --no_emp           disable calculation of empirical P-value
    --seed     -s      seed for permutations
    
    == OPTIONS FOR SELF-CONTAINED TEST ==
    --verbose  -v      print PLINK output 
    --pheno    -p      alternate phenotype files
    --out      -o      prefix for output file
    --logistic         perform logistic association test (in PLINK)
    --linear           perform logistic association test (in PLINK)
    --covar            covariate file
    --gene_based       sets are based on single genes
    --adjust           adjust pvalues of real data with lambda
    --no_plots         do not create plots within R


    === MERGE PERMUTATIONS ===
    merge multiple output files when step 2 is split in multiple processes
    --merge            merge multiple permutation files
    --out              set a prefix to find the permuted results (optional)
    
    === GET RANDOM DRAWS ===
    draw multiple random SNP or gene sets out of list with all genes and SNPs
    obligatory 
    --ndraw            number of sets of random genes/SNP sets to be drawn (= integer )
    --set              geneset file
    --pool         list with all genes and SNPs
    
    ==switch==
    --draw_ngenes           select functional set by set name
    --draw_neff_all         {set name}        select all SNPs default
    --draw_neff_genic       {set name}        only select SNPs located within genes
    --draw_neff_intergenic  {set name}        only select SNPs located outside genes
 

    obligatory with --neff option:
    --neff_calc
    --bfile
    
    mandatory: 
    --include    

    === CALCULATE EMPP OF EMPP ===
    --orig_empp        file containing the empP on the real data
    --control_empp     file containing the empP on the random draws
    --gene_set       select functional set by set name

   ===examples===:
   snp2gene annotation       jag --snp2gene --gene_loc hg18/hg18.gene.loc --snp_loc hg18/hg18.snp.loc 
   self-contained test:      jag --bfile example_data/example --set jag.set.annot --perm 10
   merge data                jag --merge test --out test
   random draws              jag --ndraw 20 --bfile example_data/example --set jag.set.annot --draw_neff_genic hsa04080 --neff_calc jag.P1.empp --pool jag.allgenes.annot
                             jag --ndraw 50 --draw_ngenes hsa04080 --set jag.set.annot --pool jag.allgenes.annot
   competitive test          jag --orig_empp jag.merged.P1.empp --control_empp jag.draw_neff_genic.P1.empp --gene_set hsa04080
    """)



class Main:
    """
    set options and runs program
    """
   
    def _step0_annotate_genes(self, out, annotate_file, up_size, down_size, gene_loc, snp_loc):
        """
        --snp2gene file with geneset(s)
        --up       size upstream region in kb
        --down     size downstream region in kb
        --gene_loc    file with gene boundaries
        --snp_loc    file with snps
        """
    
        import jag.snp2geneset as snp2geneset
        all_gene2snp, geneset2snp = snp2geneset.snp2geneset(annotate_file, up_size, down_size, gene_loc, snp_loc)
        
        set_as_string = snp2geneset.format_mapping(geneset2snp)
        all_genes_as_string = snp2geneset.format_mapping(all_gene2snp)
        
        out_sets = out.save_text_to_filename("set.annot", set_as_string)
        out_all = out.save_text_to_filename("allgenes.annot", all_genes_as_string)
   
        log.info("Saved snps2allgenes file as " + out_all)
        log.info("Saved snps2geneset file as " + out_sets)


    def _step1_run_association_analysis(self, geneset, plink, inoutput, gene_based, adjusted):
        """
        run association analysis
        """          
                             
        log.info("Running association analysis using PLINK with command:")
                
        plink.print_command = True
        from jag.association_analysis import AssociationAnalysis
        association_analysis = AssociationAnalysis(inoutput)
        files = association_analysis.run_step1(plink, geneset, gene_based, adjusted)
  
      
    def _step2_run_permutations(self, geneset, perm, seed, no_emp, plink, inoutput, gene_based): 
        """
        run permutations
        """    
        
        log.info("\nRunning permutations...")
        plink.print_command = False       
        from jag.permutation import Permutation
        permutation = Permutation(plink, geneset, inoutput)
        permutation.no_emp = no_emp
        permutation.localpermutation(perm, seed, gene_based) 
        
        if perm < 2:
            log.info("\nNo distribution plot is generated since you should run at least 2 permutations to do so.")
            inoutput.run_rproject = False
        
        
    def _step3_merge_results(self, inoutput):
        """
        merge multiple files of permutation results
        
        """
      
        import jag.mergepermutations as mergepermutations
        merger = mergepermutations.MergePermutation(inoutput)
        merger.mergeresults()
       
        log.info(_get_endtime())
        sys.exit(0)


    def _step4a_draw_genes(self, draw , draw_ngenes, ndraws):
        """
        draw random sets of genes, based on number of genes in geneset of interest
        """
    
        n = str(ndraws)
        log.info("Drawing " + n + " random gene sets...")
        import jag.drawrandom as drawrandom
        genes = drawrandom.Genes(draw)
        genes.genes(draw_ngenes, ndraws)
        log.info(_get_endtime())
        sys.exit(0)


    def _step4b_draw_snps(self, draw, neff, ndraws, draw_selection, plink, emp, out):
        """
        draw random sets of snps (based on calculated effective snps)
        """
        common.getfile_handle(emp, "(--neff_calc)", False)

        import jag.drawrandom as drawrandom
        log.info ("Drawing random SNP sets (based on nEff)...")
       
        snps = drawrandom.Snp(draw)
        snps.inregion = draw_selection

        snps.snp(neff, ndraws, emp, plink, out)
        log.info(_get_endtime())
        sys.exit(0)

    
    def _calc_emp_p_of_emp_p(self, orig, sims, set_name):
        """
        calculation of empirical pvalue for competitive test
        """
        #check if all 3 variables are set
        if sum([bool(item)for item in [orig, sims, set_name]]) == 3:
            import jag.permutationresults as permutation
            orig_header, orig_data = permutation.read_permutated_results(orig)
            try:
                row_setname = orig_data[0].index(set_name.strip())

            except ValueError:
                log.info("Could not find " + set_name + " in  the file " + orig)
                log.info(_get_terminated_time())
                sys.exit(1)

            pvalue_orig = orig_data[orig_header.index("emp_p")][row_setname]

            sims_header, sims_data = permutation.read_permutated_results(sims)
            emp_p_sims = False
            try:
                emp_p_sims = sims_data[sims_header.index("emp_p")]
            except ValueError:
                log.info("Could not find emp_p in the file" + sims)
                log.info(_get_terminated_time())
                sys.exit(1)


            log.info("Empirical P-value for competitive test for gene-set " + set_name + " = " + str(permutation.calc_emp_emp(pvalue_orig, emp_p_sims)))
            log.info(_get_endtime())

        else:
            log.info("Original empp or set is not set")
            log.info(_get_terminated_time())
            sys.exit(1)
            
    
    def extract_variables_from_command_line(self, args):
        """
        get variables from command line
        """
        
        # if no variables > print usage 
        if len(args) == 0:
            print _get_header_of_program()
            print _get_starttime()
            print "\nNo arguments are given to run JAG."
            print "Type jag --help for usage."
            print _get_terminated_time()
            sys.exit(2)
                
        else:
            try:
                opts, args = getopt.getopt(args, "p:chb:m:vg:o:s:", ["ndraw=", "draw_ngenes=", "neff_calc=",
                    "out=", "snp2gene=", "draw_neff_all=", "draw_neff_intergenic=",
                    "draw_neff_genic=", "exclude",
                    "pool=", "snp_loc=", "gene_loc=", "down=", "up=",
                    "include", "no_emp", "gene_based", "no_plots", "help",
                    "merge=",
                    "verbose", "runplink=",
                    "bfile=", "pheno=", "set=", "perm=", "seed=",
                    "orig_empp=", "control_empp=", "gene_set=",
                    "covar=", "logistic", "linear", "adjust"])
            
            except getopt.GetoptError, err:
                print _get_header_of_program()
                print _get_starttime()
                print ""
                print str(err)
                print "Type jag --help for usage."
                print _get_terminated_time()
                sys.exit(2)
                
        return opts


    def enable_logging_with_prefix(self, inoutput):
        """
        #log with prefix
        """
        log.basicConfig(level=log.INFO, stream=sys.stdout, #format='%(levelname)s\t%(asctime)s   %(filename)s:%(lineno)d\t%(message)s',
            format='%(message)s', datefmt='%H:%M:%S')
        file_handler = log.FileHandler(inoutput.out + 'log', mode='w')
        log.root.addHandler(file_handler)
        

    def __init__(self, args):
                
        opts = self.extract_variables_from_command_line(args)
        
        from jag.plink import Plink
        plink = Plink()
                       
        geneset = None
        perm = None
        seed = None

        annotate_file = False
        up = float(0)
        down = float(0)
        gene_loc_file = None
        snp_loc = None

        no_emp = False
        create_plots = True
        draw_ngenes = False
        empirical_p_filename = False

        ndraws = 0
        draw_selection = None
        exclude_group = True
        complete_gene_snp_mapping = None
        verbose = False
        gene_based = False
       
        merge = False
        adjusted = False
        prefix = None
        
        orig = False
        sims = False
        set_name = False

        out_prefix = "jag"
        
        # open log file
        for identifier, assigned_value in opts:
            if identifier in ("-o", "--out"):
                out_prefix = assigned_value
                
            elif identifier in ("-s", "--set"):
                match_draw = re.match("\w+(\.draws_n\w*)\.set.annot$", assigned_value)
                if match_draw is not None:
                    out_prefix = out_prefix + match_draw.group(1)
                 
            elif identifier in ("--gene_based"):
                gene_based = True
                out_prefix = out_prefix + ".gene_based"
              
        from jag.file_fetch_and_write import InAndOut
        inoutput = InAndOut()
        inoutput.set_outfile(out_prefix)
        
        self.enable_logging_with_prefix(inoutput)        
        log.info(_get_header_of_program())
        log.info("Save logfile as [" + inoutput.out + "log]\n")
        log.info(_get_starttime())
        
        log.info("\nUsed options:")
        for o, a in opts:
            if o not in ("-h", "--help"):
                log.info("\t" + o + " " + a)
            else:
                log.info("\t" + o + " " + a)
                log.info("\nPrinting help documentation...")
                log.info(_usage())
                print(_get_terminated_time())
                sys.exit(2)
     
        if len(opts) == 0:
            _usage()
            
        else:
            
            for identifier, assigned_value in opts:

                if identifier in ("-o", "--out"):
                    #this stuff is printed to catch the out for the location of the log file
                    out_prefix = assigned_value
                                       
                elif identifier in ("-s", "--set"):
                    assert common.getfile_handle(assigned_value,"(--set)", verbose)
                    group = assigned_value
                    plink.group = assigned_value
                    
                elif identifier in ("-m", "--perm"):
                    perm = int(assigned_value)
                    
                elif identifier in ("-v", "--verbose"):
                    plink.verbose = True
                    verbose = True
                    
                elif identifier in ("--snp_loc"):
                    snp_loc = assigned_value
                    
                    assert common.getfile_handle(assigned_value, \
                                                 "(--snp_loc)", verbose)

                elif identifier in ("--control_empp"):
                    assert common.getfile_handle(assigned_value, \
                                                  "(--control_empp)", verbose)
                    sims = assigned_value
                    
                elif identifier in ("--orig_empp"):
                    assert common.getfile_handle(assigned_value, \
                                                  "(--orig_empp)", verbose)
                    orig = assigned_value
                    
                elif identifier in ("--gene_set"):
                    set_name = assigned_value
                    
                elif identifier in ("--no_emp"):
                    no_emp = True
                    create_plots = False
                    
                elif identifier in ("--no_plots"):
                    create_plots = False
                    
                elif identifier in ("--ndraw"):
                    ndraws = int(assigned_value)
                    
                elif identifier in ("--draw_ngenes"):
                    draw_ngenes = assigned_value
                    
                elif identifier in ("--snp2gene"):
                    assert common.getfile_handle(assigned_value, \
                                                 "(--snp2gene)", verbose)
                    annotate_file = assigned_value
                    
                elif identifier in ("--up"):
                    try:
                        up = float(assigned_value)
                        
                    except ValueError:
                        log.info("Value after --up parameter should be a number")
                        sys.exit(1)
                                       
                elif identifier in ("--down"):
                    try:
                        down = float(assigned_value)
                        
                    except ValueError:
                        log.info("Value after --down parameter should be a number")
                        sys.exit(1)
                
                elif identifier in ("--gene_loc"):
                    assert common.getfile_handle(assigned_value, "Cannot find file with gene boundaries (set by --gene_loc)", verbose)
                    gene_loc_file = assigned_value
                 
                elif identifier in ("--pool"):
                    complete_gene_snp_mapping = assigned_value

                elif identifier in ("--draw_neff_genic"):
                    draw_selection = "genic"
                    set_name = assigned_value
                    
                elif identifier in ("--draw_neff_intergenic"):
                    draw_selection = "intergenic"
                    set_name = assigned_value
                    
                elif identifier in ("--draw_neff_all"):
                    draw_selection = "all"
                    set_name = assigned_value
                    
                elif identifier in ("--covar"):
                    assert common.getfile_handle(assigned_value, \
                         "(--covar;", verbose)
                    plink.covar_file = assigned_value
                                
                elif identifier in ("--linear"):
                    plink.switches += "--linear "
                                        
                elif identifier in ("--logistic"):
                    plink.switches += "--logistic "
                    
                elif identifier in ("--adjust"):
                    plink.switches += "--adjust "
                    adjusted = True

                elif identifier in ("--exclude"):
                    exclude_group = True
                    
                elif identifier in ("--include"):
                    exclude_group = False

                elif identifier in ("--neff_calc"):
                    assert common.getfile_handle(assigned_value, "(--neff_calc)", verbose)
                    empirical_p_filename = assigned_value
                    
                elif identifier in ("-b", "--bfile"):
                    plink.bfile = common.check_bim_bed_fam(assigned_value)
                    
                elif identifier in ("-p", "--pheno"):
                    assert common.getfile_handle(assigned_value, \
                     "(-p or --pheno)", verbose)
                    plink.pheno_file = assigned_value
                    
                elif identifier in ("--seed"):
                    seed = assigned_value
                    
                elif identifier in ("--gene_based"):
                    gene_based = True

                elif identifier in ("--merge"):
                    prefix = assigned_value
                    inoutput = InAndOut()
                    inoutput.set_outfile(prefix)
                    merge = True
                                    
                else:
                    assert False, "unhandled option"
                    print _get_terminated_time()
                    sys.exit(2)
                    
            log.info("")        
                  
            inoutput.run_rproject = create_plots
                       
            if merge:
                self._step3_merge_results(inoutput)

            if orig or sims:
                self._calc_emp_p_of_emp_p(orig , sims , set_name)
                sys.exit(0)

            if annotate_file and gene_loc_file and snp_loc:
                self._step0_annotate_genes(inoutput, annotate_file, up, down, \
                                            gene_loc_file, snp_loc)

            elif draw_ngenes or draw_selection:
                common.getfile_handle(complete_gene_snp_mapping, "--pool is not set correctly.", verbose)
                common.getfile_handle(group, "--set is not set.", verbose)

                import jag.drawrandom as drawrandom
                draw = drawrandom.DrawRandom(complete_gene_snp_mapping, group, inoutput)

                draw.exclude = exclude_group

                if seed:
                    draw.setseed(seed)

                if draw_ngenes:
                    self._step4a_draw_genes(draw, draw_ngenes, ndraws)
                
                if draw_selection:
                    self._step4b_draw_snps(draw, set_name, ndraws, draw_selection, \
                                           plink, empirical_p_filename, out_prefix)

            elif(perm > 0 and plink.bfile and group):
                #Run Step 1 + Step 2
                self._step1_run_association_analysis(group, plink, inoutput, gene_based, adjusted)
                self._step2_run_permutations(group, perm, seed, no_emp, plink, inoutput, gene_based)

            elif(perm == 0):
                # Run only association analysis
                self._step1_run_association_analysis(group, plink, inoutput, gene_based, adjusted)
                log.info("\nNo permutations will be proceeded since number of permutations is zero.")
                
            else:
                log.info("You are using an invalid combination of parameters.\nCheck the help file for the correct combination.")

        log.info(_get_endtime())
            
if __name__ == '__main__':
    Main(sys.argv[1:])
