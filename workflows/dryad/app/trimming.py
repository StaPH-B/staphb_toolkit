#author: Kelsey Florek
#email: kelsey.florek@slh.wisc.edu
#quality trimming function

import os, sys

#local libs
from workflows.dryad.app.lib import checkexists
import core.sb_programs

def q_trim(reads,jobs,cpu,outdir,tracker):
    minlength = 100
    windowsize = 4
    qscore = 30
    logfile = os.path.join(outdir,'qtrim.log')

    cmds = []
    read_path = ''
    for read_pair in reads:
        #main command
        main_cmd = 'java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads {0}'.format(cpu)

        sid = os.path.basename(read_pair[0]).split('_')[0]

        args = ' {read1} {read2} -baseout /output/{sid}.fastq.gz SLIDINGWINDOW:{windowsize}:{qscore} MINLEN:{minlength}'.format(minlength=minlength,windowsize=windowsize,qscore=qscore,read1=os.path.basename(read_pair[0]),read2=os.path.basename(read_pair[1]),sid=sid)

        cmds.append(main_cmd + args)

        if read_path == '':
            read_path = os.path.dirname(os.path.abspath(read_pair[1]))
        elif read_path != os.path.dirname(os.path.abspath(read_pair[1])):
            print("Reads cannot be in multiple locations. Exiting.")
            sys.exit()
        else:
            pass
    checkexists(os.path.join(outdir,"trimmed"))

    #start multiprocessing
    pool = mp.Pool(processes=jobs)
    print("Begining quality trimming of reads:\n Number of Jobs: {0}\n CPUs/Job: {1}".format(jobs,cpu))
    #denote logs
    with open(logfile,'a') as outlog:
        outlog.write('***********\n')
        outlog.write('Trimmomatic\n')
        #begin multiprocessing
        results = pool.starmap_async(cd.call,[['staphb/trimmomatic:0.39',cmd,'/data',{read_path:"/data",os.path.join(outdir,'trimmed'):"/output"}] for cmd in cmds])
        stdouts = results.get()
        for stdout in stdouts:
            outlog.write('-----------\n')
            outlog.write(stdout)
        #denote end of logs
        outlog.write('***********\n')

    #remove unpaired reads
    for root,dirs,files in os.walk(os.path.join(outdir,'trimmed')):
        for file in files:
            if "U.fastq.gz" in file:
                os.remove(os.path.join(root,file))

    print("Finished Quality Trimming Reads")

    #update status
    tracker.update_status_done('trimmed')
