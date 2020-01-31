import os,sys
import psutil
import shutil
from shutil import which
import json
import datetime
import getpass
import signal
import csv

from staphb_toolkit.core import fileparser
from staphb_toolkit.core.sb_programs import Run
from staphb_toolkit.core.sb_programs import Run_multi
from staphb_toolkit.lib import sb_mash_species

from workflows.spriggan.app.lib import StatusTracker as ST
from workflows.spriggan.app.lib import getfiles,checkexists,cpu_count,run_and_log

def species(read_path,config_path,jobs,cpu_job,spriggan_out):
    # set up Mash species
    mash_obj = sb_mash_species.MashSpecies(path=read_path, output_dir=spriggan_out, configuration=config_path)

    # if we don't have mash species completed run it, otherwise parse the file and get the results
    mash_results = os.path.join(*[spriggan_out, 'mash_output', 'mash_species.csv'])
    if not os.path.isfile(mash_results):
        print(f'Making taxonomic predictions')
        mash_species = mash_obj.run()

    else:
        mash_species = {}
        with open(mash_results, 'r') as csvin:
            reader = csv.reader(csvin, delimiter=',')
            for row in reader:
                mash_species[row[0]] = row[1]

# trimming function
def q_trim(read_path,read_dict,spriggan_config,jobs,cpu_job,spriggan_out):
    # trimming output path
    trimmed_path = os.path.join(spriggan_out,'trimmed')
    input_path = os.path.join(spriggan_out,'input_reads')
    checkexists(trimmed_path)

    # create path to log file
    logfile = os.path.join(spriggan_out,'qtrim.log')

    # set up params, docker mounting, image and tag
    trimmomatic_configuration = spriggan_config['parameters']['trimmomatic']

    trimmomatic_params = trimmomatic_configuration['params']
    windowsize = trimmomatic_params['windowsize']
    qscore = trimmomatic_params['qscore']
    minlength = trimmomatic_params['minlength']

    trimmomatic_mounting = {read_path:'/data',trimmed_path:'/output'}
    image = trimmomatic_configuration['image']
    tag = trimmomatic_configuration['tag']

    # set up Trimmomatic command list
    cmds = []

    for id in read_dict:
        # sample IDs
        read1 = os.path.basename(read_dict[id].fwd)
        read2 = os.path.basename(read_dict[id].rev)
        sid = read1.split('_')[0]

        # set up cmds
        main_cmd = f'java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads {cpu_job}'
        args = f' {sid}/{read1} {sid}/{read2} -baseout /output/{sid}.fastq.gz SLIDINGWINDOW:{windowsize}:{qscore} MINLEN:{minlength}'
        cmds.append(main_cmd + args)

    print(f'Beginning quality trimming reads:\n Number of Jobs: {jobs}\n CPUs/Job: {cpu_job}')

    # create run object and run commands
    process = 'Trimming with Trimmomatic'
    runs = Run_multi(command_list=cmds,path=trimmomatic_mounting,image=image,tag=tag)
    run_and_log(process,logfile,runs,jobs)

    # remove unpaired reads
    for root,dirs,files in os.walk(trimmed_path):
        for file in files:
            if 'U.fastq.gz' in file:
                os.remove(os.path.join(root,file))

    print('Finished quality trimming reads')

    # update config file
    spriggan_config['file_io']['output_files']['log_files']['qtrim'] = os.path.abspath(os.path.join(spriggan_out,'qtrim.log'))

# assembly function
def assemble_reads(read_dict,spriggan_config,jobs,cpu_job,spriggan_out):
    # get trimmed reads
    trimmed_path = os.path.join(spriggan_out,'trimmed')
    trimmed_dict = read_dict

    for id in trimmed_dict:
        trimmed_dict[id].path = trimmed_path
        trimmed_dict[id].fwd = os.path.basename(read_dict[id].fwd).replace("_1.fastq.gz", "_1P.fastq.gz")
        trimmed_dict[id].rev = os.path.basename(read_dict[id].rev).replace("_2.fastq.gz", "_2P.fastq.gz")

    # update config file output
    if spriggan_config['file_io']['output_files']['trimmed_fastqs'] == None:
        spriggan_config['file_io']['output_files']['trimmed_fastqs'] = {}

    for id in trimmed_dict:
        spriggan_config['file_io']['output_files']['trimmed_fastqs'][trimmed_dict[id].id] = [os.path.join(trimmed_path,trimmed_dict[id].fwd),os.path.join(trimmed_path,trimmed_dict[id].rev)]

    # assembly output path
    assemblies_path = os.path.join(spriggan_out,'assemblies')
    checkexists(assemblies_path)

    # create log file
    logfile = os.path.join(spriggan_out,'assembly.log')

    # determine free ram
    free_ram = int(psutil.virtual_memory()[1]/1000000000)
    ram_job = int(free_ram / jobs)

    # set up params, Docker mounting, image and tag
    shovill_configuration = spriggan_config['parameters']['shovill']
    shovill_params = shovill_configuration['params']

    shovill_mounting = {trimmed_path:'/data',assemblies_path:'/output'}
    image = shovill_configuration['image']
    tag = shovill_configuration['tag']

    # set up Shovill command list
    cmds = []

    for id in trimmed_dict:
        # sample IDs
        read1 = os.path.basename(trimmed_dict[id].fwd)
        read2 = os.path.basename(trimmed_dict[id].rev)
        sid = read1.split('_')[0]

        # set up cmds
        cmds.append(f'bash -c \"shovill --R1 /data/{read1} --R2 /data/{read2} --force --outdir /output/{sid} --cpus {cpu_job} --ram {ram_job} && mv /output/{sid}/contigs.fa /output/{sid}/{sid}.fa\"')

    print(f'Beginning assembly of reads:\n Number of Jobs: {jobs}\n CPUs/Job: {cpu_job}')

    # create run object and run commands
    process = 'Assembly with Shovill'
    runs = Run_multi(command_list=cmds,path=shovill_mounting,image=image,tag=tag)
    run_and_log(process,logfile,runs,jobs)

    print('Finished assembling reads')

    # update config file
    if spriggan_config['file_io']['output_files']['assemblies'] == None:
        spriggan_config['file_io']['output_files']['assemblies'] = {}

    for id in trimmed_dict:
        sid = os.path.basename(trimmed_dict[id].fwd.split('_')[0])
        spriggan_config['file_io']['output_files']['assemblies'][sid] = os.path.join(assemblies_path,f'{sid}/{sid}.fa')

    spriggan_config['file_io']['output_files']['log_files']['assembly'] = os.path.join(spriggan_out,'assembly.log')

def mlst(spriggan_config,jobs,cpu_job,spriggan_out):
    # get assemblies
    assemblies_path = os.path.join(spriggan_out,'assemblies')
    fastas = list(getfiles(assemblies_path))[0]

    # mlst output path
    mlst_path = os.path.join(spriggan_out,'mlst')
    checkexists(mlst_path)

    # set up params, docker mounting, image and tag
    mlst_configuration = spriggan_config['parameters']['mlst']
    mlst_params = mlst_configuration['params']

    mlst_mounting = {assemblies_path:'/data',mlst_path:'/output'}
    image = mlst_configuration['image']
    tag = mlst_configuration['tag']

    # set up string of Shovill results, \s delimited
    fasta_files = ''
    for path in fastas:
        assembly_file = os.path.basename(path)
        assembly_dir = os.path.basename(os.path.dirname(path))
        fasta_files = fasta_files + f'/data/{assembly_dir}/{assembly_file} '

    # main command
    cmd = f'bash -c \"mlst --nopath --quiet ' + fasta_files + '> /output/mlst.tab\"'

    print('Predicting MLST types')

    # set up run and run process
    mlst_run = Run(command=cmd,path=mlst_mounting,image=image,tag=tag)
    mlst_run.run()

    print('Finished predicting MLST types')

    # copy output to main dir
    shutil.copyfile(os.path.join(mlst_path,'mlst.tab'),os.path.join(spriggan_out,'mlst.tab'))

    # update config file
    spriggan_config['file_io']['output_files']['mlst'] = os.path.join(spriggan_out,'mlst.tab')

# AR prediction function
def abricate(spriggan_config,jobs,cpu_job,spriggan_out):
    # get assemblies
    assemblies_path = os.path.join(spriggan_out,'assemblies')
    fastas = list(getfiles(assemblies_path))[0]

    # Abricate output path
    abricate_path = os.path.join(spriggan_out,'abricate')
    checkexists(abricate_path)

    # set up params, docker mounting, image and tag
    abricate_configuration = spriggan_config['parameters']['abricate']
    abricate_mounting = {assemblies_path:'/data',abricate_path:'/output'}
    image = abricate_configuration['image']
    tag = abricate_configuration['tag']

    # run Abricate on assemblies
    print('Predicting AR')

    # set up commands for Abricate
    cmds = []
    for path in fastas:
        # main command
        assembly_file = os.path.basename(path)
        assembly_dir = os.path.basename(os.path.dirname(path))
        sid = assembly_file.split('.')[0]
        handle = sid + '.tab'
        cmds.append(f'bash -c \"abricate --db ncbi --quiet /data/{assembly_dir}/{assembly_file} > /output/{handle}\"')

    # set up run and job number
    abricate_run = Run_multi(command_list=cmds,path=abricate_mounting,image=image,tag=tag)
    njobs = jobs
    abricate_run.run(jobs=njobs)

    print('Finished predicting AR')

    # update config file output
    if spriggan_config['file_io']['output_files']['abricate'] == None:
        spriggan_config['file_io']['output_files']['abricate'] = {}

    for path in fastas:
        sid = os.path.basename(path).split('.')[0]
        spriggan_config['file_io']['output_files']['abricate'][sid] = os.path.join(abricate_path,sid + '.tab')

    # summarize Abricate results
    print('Summarizing AR prediction')

    abricate_mounting = {abricate_path:'/data',spriggan_out:'/output'}

    # get Abricate results
    tsvs = list(getfiles(abricate_path))[0]

    # set up string of Abricate results, \s delimited
    tsv_files = ''
    for path in tsvs:
        abricate_file = os.path.basename(path)
        tsv_files = tsv_files + f'/data/{abricate_file} '

    # main command
    cmd = f'bash -c \"abricate --summary ' + tsv_files + '> /output/abricate_summary.tab\"'

    # set up run
    summary_run = Run(command=cmd,path=abricate_mounting,image=image,tag=tag)
    summary_run.run()

    print('Finished summarizing AR prediction')

    # update config file output
    spriggan_config['file_io']['output_files']['abricate_summary'] = os.path.abspath(os.path.join(spriggan_out,'abricate_summary.tab'))

# assembly QC function
def quast(reference,spriggan_config,jobs,cpu_job,spriggan_out):
    # get assemblies
    assemblies_path = os.path.join(spriggan_out,'assemblies')
    fastas = list(getfiles(assemblies_path))[0]

    # path to reference genome
    ref_path = os.path.dirname(os.path.abspath(reference))
    ref_genome = os.path.basename(os.path.abspath(reference))

    # quast output path
    quast_path = os.path.join(spriggan_out,'quast')
    checkexists(quast_path)

    # create path to log file
    logfile = os.path.join(spriggan_out,'quast.log')

    # set up params, docker mounting, image and tag
    quast_configuration = spriggan_config['parameters']['quast']
    quast_params = quast_configuration['params']

    quast_mounting = {assemblies_path:'/data',ref_path:'/reference',quast_path:'/output'}
    image = quast_configuration['image']
    tag = quast_configuration['tag']

    # add reference genome to config file
    spriggan_config['parameters']['quast']['params']['reference'] = ref_genome

    # set up string of Shovill results, \s delimited
    fasta_files = ''
    for path in fastas:
        assembly_file = os.path.basename(path)
        assembly_dir = os.path.basename(os.path.dirname(path))
        fasta_files = fasta_files + f'/data/{assembly_dir}/{assembly_file} '

    # main command
    cmd = f'quast.py --silent -t {cpu_job} -o /output/ -r /reference/{ref_genome} '+ fasta_files

    # set up run and run process
    quast_run = Run(command=cmd,path=quast_mounting,image=image,tag=tag)
    quast_run.run()

    print('Finished running Quast')

    # copy output to main dir
    shutil.copyfile(os.path.join(quast_path,'transposed_report.tsv'),os.path.join(spriggan_out,'quast_report.tsv'))
    shutil.copyfile(os.path.join(quast_path,'quast.log'),os.path.join(spriggan_out,'quast.log'))

    # update config file
    spriggan_config['file_io']['output_files']['quast_report'] = os.path.join(spriggan_out,'quast_report.tsv')
    spriggan_config['file_io']['output_files']['log_files']['quast'] = os.path.join(spriggan_out,'quast.log')

def spriggan(read_path,reference,threads,outdir,spriggan_config):
    # set cwd as output dir if outdir is empty
    try:
        outdir = os.path.abspath(outdir)
    except (AttributeError,TypeError) as err:
        outdir = os.getcwd()

    # get num of jobs and number of cpus per job
    jobs,cpu_job = cpu_count(threads)

    # initialize tracker
    pipeline = 'spriggan'
    tracker = ST()
    spriggan_out = tracker.initialize(outdir,pipeline)
    project = os.path.basename(spriggan_out)

    # get fastqs
    fastq_files = fileparser.ProcessFastqs(read_path)

    # set read_path equal to outdir after reads have been copied there
    fastq_files.link_reads(output_dir=spriggan_out)

    # path to untrimmed reads
    read_path = os.path.join(spriggan_out, 'input_reads')
    read_dict = fastq_files.id_dict()

    # get config file
    if spriggan_config is not None:
        config_path = os.path.abspath(spriggan_config)
    else:
        # use default
        config_path = os.path.join(os.path.abspath(os.path.dirname(os.path.realpath(__file__))),'spriggan_config.json')

    # load config file
    with open(config_path) as config_file:
        spriggan_config = json.load(config_file)

    # update configuration object with fastqs, reference, run ID, user and date/time
    spriggan_config['file_io']['input_files']['reference'] = os.path.abspath(reference)
    spriggan_config['execution']['run_id'] = project
    spriggan_config['execution']['user'] = getpass.getuser()
    spriggan_config['execution']['date_time'] = datetime.datetime.today().strftime('%Y-%m-%d')

    if spriggan_config['file_io']['input_files']['fastqs'] == None:
        spriggan_config['file_io']['input_files']['fastqs'] = {}

    for fastqs in fastq_files.reads:
        spriggan_config['file_io']['input_files']['fastqs'][fastq_files.reads[fastqs].id] = [fastq_files.reads[fastqs].fwd,fastq_files.reads[fastqs].rev]

    # begin wokflow
    if not tracker.check_status('species'):
        print('Predicting species using Mash')
        species(read_path,config_path,jobs,cpu_job,spriggan_out)
        tracker.update_status_done('species')

    if not tracker.check_status('trimmed'):
        print('Trimming reads using Trimmomatic')
        q_trim(read_path,read_dict,spriggan_config,jobs,cpu_job,spriggan_out)
        tracker.update_status_done('trimmed')

    if not tracker.check_status('assemble'):
        print('Assembling reads using Shovill')
        assemble_reads(read_dict,spriggan_config,jobs,cpu_job,spriggan_out)
        tracker.update_status_done('assemble')

    if not tracker.check_status('mlst'):
        print('Predicting MLST types using MLST')
        mlst(spriggan_config,jobs,cpu_job,spriggan_out)
        tracker.update_status_done('mlst')

    if not tracker.check_status('abricate'):
        print('IDing AR genes using Abricate')
        abricate(spriggan_config,jobs,cpu_job,spriggan_out)
        tracker.update_status_done('abricate')

    if not tracker.check_status('quast'):
        print('Evaluating assembly quality using Quast')
        quast(reference,spriggan_config,jobs,cpu_job,spriggan_out)
        tracker.update_status_done('quast')

    # write json file to output dir
    spriggan_config_file = os.path.join(spriggan_out, project+'_config.json')
    with open(spriggan_config_file, 'w') as jsonOut:
        jsonOut.write(json.dumps(spriggan_config, indent=3))
