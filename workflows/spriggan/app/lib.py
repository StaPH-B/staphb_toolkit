import os,sys
import time
import datetime
from staphb_toolkit.core.sb_programs import Run_multi
from staphb_toolkit.core.sb_programs import Run

# returns (num jobs,num cpus per job)
def cpu_count(num):
    if num <= 1:
        return 1,1
    elif num <= 5:
        return 1,num
    else:
        results = []
        for n in range(2,num):
            if num % n == 0:
                results.append(n)
        if len(results) < 2:
            return cpu_count(num-1)
        index = int(len(results)/2)
        return results[index],int(num/results[index])

def checkexists(path):
    path = os.path.abspath(path)
    if not os.path.isdir(path):
        os.mkdir(path)
        return False
    else:
        return True

# object for tracking the progress of the pipeline
class StatusTracker:
    # status code defined
    # 0 = created directory, ready to start
    # 1 = species
    # 2 = trimmed
    # 3 = assemble
    # 4 = mlst
    # 5 = abricate
    # 6 = quast
    def __init__(self):
        self.status_code = {'start':'0','species':'0','trimmed':'0','assemble':'0','mlst':'0','abricate':'0','quast':'0'}
        self.status_file_path = ''

    def reset_code(self):
        self.status_code = {'start':'0','species':'0','trimmed':'0','assemble':'0','mlst':'0','abricate':'0','quast':'0'}

    # check if we have completed the requested pipeline, internal function
    def checkComplete(self,pipeline):
        if self.status_code['quast'] == '1':
            return True
        else:
            return False

    def write_status_code(self):
        code_string = ''
        for key in self.status_code:
            code_string += self.status_code[key]
        with open(self.status_file_path,'w') as outstat:
            outstat.write(code_string)

    def read_status_code(self):
        with open(self.status_file_path,'r') as instat:
            code_string = instat.readline()
        try:
            c = 0
            for key in self.status_code:
                self.status_code[key] = code_string[c]
                c += 1
        except:
            print("Saved status code is longer than possible")
            sys.exit(2)

    # update the status on a finished process
    def update_status_done(self,process):
        self.status_code[process] = '1'
        self.write_status_code()

    # check the status of completion on a process
    def check_status(self,process):
        self.read_status_code()
        if self.status_code[process] == '1':
            return True
        else:
            return False

    # initalize the tracker object
    def initialize(self,path,pipeline):
        # check for a previous run
        for root,dirs,files in os.walk(path):
            for dir in dirs:
                if "spriggan-" in dir:
                    spriggan_path = os.path.join(root,dir)
                    self.status_file_path = os.path.join(spriggan_path,'tracker')
                    try:
                        self.read_status_code()
                    except:
                        print("Cannot read tracker file please delete:")
                        print(spriggan_path)
                        print("and try again.")
                        sys.exit(2)
                    if not self.checkComplete(pipeline):
                        print("There is a previous unfinished run, do you wish to continue?")
                        print(spriggan_path)
                        continue_run = input("Y/N? ")
                        if continue_run == 'n' or continue_run == 'N':
                            pass
                        elif continue_run == 'y' or continue_run == 'Y':
                            self.update_status_done(pipeline)
                            return spriggan_path
                        else:
                            print("Not a Y/N, exiting!")
                            sys.exit(2)

        # create output dir and new tracker
        self.reset_code()
        time = datetime.datetime.now()
        str_time = "{0}{1}{2}{3}{4}".format(time.year,time.month,time.day,time.hour,time.minute)
        if checkexists(os.path.join(path,"spriggan-"+str_time)):
            print(os.path.join(path,"spriggan-"+str_time)+" already exists... please try again in 1 min.")
            sys.exit(2)
        spriggan_path = os.path.join(path,"spriggan-"+str_time)
        self.status_file_path = os.path.join(spriggan_path,'tracker')
        self.update_status_done('start')
        self.update_status_done(pipeline)
        return spriggan_path


def getfiles(path):
    # scan path and look for files
    fastq_files = []
    fasta_files = []
    tsv_files = []

    for root,dirs,files in os.walk(path):
        for file in files:
            if ".fastq.gz" in file:
                fastq_files.append(os.path.join(root,file))
            if ".fa" in file and "spades" not in file:
                fasta_files.append(os.path.join(root,file))
            if ".tab" in file:
                tsv_files.append(os.path.join(root,file))

    if len(fastq_files) > 0:
        fastq_files.sort()
        if len(fastq_files) % 2 != 0:
            print('There is an uneven number of read pairs in {0}. Exiting.'.format(path))
            sys.exit()
        paired_reads = []
        [paired_reads.append([x,y]) for x,y in zip(fastq_files[0::2],fastq_files[1::2])]
        yield paired_reads

    if len(fasta_files) > 0:
        yield fasta_files

    if len(tsv_files) > 0:
        yield tsv_files

def run_and_log(process, logfile, runs, njobs):
    # save stdout
    original = sys.stdout

    # set stdout to log file and run process
    sys.stdout = open(logfile, 'a')
    print('***********\n')
    print('-----------\n')
    print(process+"\n")
    runs.run(jobs=njobs)
    print('-----------\n')
    print('***********\n')

    # reset stdout
    sys.stdout = original
