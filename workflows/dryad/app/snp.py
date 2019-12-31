#author: Kelsey Florek
#email: kelsey.florek@slh.wisc.edu
#snp pipeline


import os,sys
import psutil
import shutil

from workflows.dryad.app.lib import getfiles,checkexists
import core.sb_programs


#CFSAN SNP pipeline
def cfsan_snp(outdir,reference,cpus):
    print("Running the CFSAN SNP Pipeline")
    reference_name = os.path.basename(reference)
    shutil.copyfile(reference,os.path.join(outdir,reference_name))
    logfile = os.path.join(outdir,'cfsan_snp.log')
    command = "sh -c 'cfsan_snp_pipeline run {0} -o /data/cfsan -s /data/snp_reads'".format(reference_name)
    #denote logs
    with open(logfile,'a') as outlog:
        outlog.write('***********\n')
        outlog.write('CFSAN SNP\n')
        stdout = cd.call('staphb/cfsan-snp-pipeline:2.0.2',command,'/data',{outdir:"/data"},cpu_set=cpus,sig_default=False)
        outlog.write('-----------\n')
        outlog.write(stdout)
        #denote end of logs
        outlog.write('***********\n')
    print("Finished CFSAN SNP Pipeline")

#create tree
def build_tree(outdir,model='GTR+G'):
    logfile = os.path.join(outdir,'tree.log')
    input_path = os.path.join(outdir,"cfsan")

    #remove tree dir if it exists
    shutil.rmtree(os.path.join(outdir,'snp_tree'),ignore_errors=True)
    #create path
    cg_path = os.path.join(outdir,"snp_tree")
    checkexists(cg_path)
    shutil.copyfile(os.path.join(input_path,'snpma.fasta'),os.path.join(cg_path,'snpma.fasta'))

    #count the number of isolates
    with open(os.path.join(os.path.join(cg_path,'snpma.fasta')),'r') as infasta:
        isolate_counter = 0
        for line in infasta.readlines():
            if line[0] == '>':
                isolate_counter += 1

    #iqtree command
    if isolate_counter > 3:
        command = "sh -c 'iqtree -s snpma.fasta -m {0} -bb 1000 '".format(model)
        tree_file = 'snpma.fasta.contree'
    else:
        command = "sh -c 'iqtree -s snpma.fasta -m {0} '".format(model)
        tree_file = 'snpma.fasta.treefile'

    print("Building the Phylogeny")

    #denote logs
    with open(logfile,'a') as outlog:
        outlog.write('***********\n')
        outlog.write('Building Tree\n')
        stdout = cd.call('staphb/iqtree:1.6.7',command,'/data',{cg_path:"/data"},sig_default=False)
        outlog.write('-----------\n')
        outlog.write(stdout)
        #denote end of logs
        outlog.write('***********\n')
    print("Finished Bulding Phylogeny")
    return tree_file


def snp(jobs,cpu_job,outdir,reference,tracker):
    #create sym links for read pairs in subfolders for cfsan snp pipeline
    if not tracker.check_status('initalize'):
        input_path = os.path.join(outdir,'trimmed')
        checkexists(os.path.join(outdir,'snp_reads'))
        fastqs = list(getfiles(input_path))[0]
        for read_pair in fastqs:
            sid = os.path.basename(read_pair[0]).split('_')[0]
            read_path = os.path.join(*[outdir,'snp_reads',sid])
            checkexists(read_path)
            os.link(read_pair[0],os.path.join(read_path,os.path.basename(read_pair[0])))
            os.link(read_pair[1],os.path.join(read_path,os.path.basename(read_pair[1])))
        tracker.update_status_done('initalize')

    if not tracker.check_status('cfsan'):
        total_cpus = jobs * cpu_job
        cfsan_snp(outdir,reference,total_cpus)
        dist_matrix_path = os.path.join(*[outdir,'cfsan','snp_distance_matrix.tsv'])
        shutil.copyfile(dist_matrix_path,os.path.join(outdir,'snp_distance_matrix.tsv'))
        tracker.update_status_done('cfsan')

    if not tracker.check_status('snp_tree'):
        tree_file = build_tree(outdir)
        in_path = os.path.join(*[outdir,'snp_tree',tree_file])
        out_path = os.path.join(outdir,'snp_tree.tree')
        shutil.copyfile(in_path,out_path)
        tracker.update_status_done('snp_tree')
