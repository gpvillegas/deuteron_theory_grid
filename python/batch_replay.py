import sys
import os

def check_dir(this_dir):
    if os.path.isdir(this_dir):
        print(this_dir, " exists, will use it ")
    else:
        # create output directory (if necessary)
        try:
            os.mkdir(this_dir)
        except Exception as err:
            print(f"problem creating directory {this_dir}: ", err)
            sys.exit()

def write_header(f,job_name,email=True):
    if email:
        f.write('#!/bin/bash\n'+\
                '#SBATCH --partition=production\n'+\
                '#SBATCH --account=hallc\n'+\
                f'#SBATCH --job-name={job_name}\n'+\
                '#SBATCH --output=/farm_out/%u/%x-%j-%N.out\n'+\
                '#SBATCH --error=/farm_out/%u/%x-%j-%N.err\n'+\
                '#SBATCH --ntasks=1\n'+\
                '#SBATCH --cpus-per-task=1\n'+\
                '#SBATCH -N1\n'+\
                '#SBATCH --mem-per-cpu=3000\n'+\
                '#SBATCH --mail-user=gvill\n'+\
                '#SBATCH --mail-type=ALL\n'+\
                '#SBATCH --time=24:00:00\n')
    else:
                f.write('#!/bin/bash\n'+\
                '#SBATCH --partition=production\n'+\
                '#SBATCH --account=hallc\n'+\
                f'#SBATCH --job-name={job_name}\n'+\
                '#SBATCH --output=/farm_out/%u/%x-%j-%N.out\n'+\
                '#SBATCH --error=/farm_out/%u/%x-%j-%N.err\n'+\
                '#SBATCH --ntasks=1\n'+\
                '#SBATCH --cpus-per-task=1\n'+\
                '#SBATCH -N1\n'+\
                '#SBATCH --mem-per-cpu=3000\n'+\
                '#SBATCH --time=24:00:00\n')            

sbatch_DIR = './current/'
check_dir(sbatch_DIR)

replay_DIR = '../deut_offline_replay/'
replay_script = 'replay_deut_prod.sh'
nevents = '-1'

with open('runs_to_replay.txt','r') as readf:
   runs_to_replay =[]
   for line in readf:
        runs_to_replay.append(line.strip())
#print (runs_to_replay)

for i in range(len(runs_to_replay)):
    replay_cmd = replay_DIR + replay_script + ' ' + f'{runs_to_replay[i]}' + ' ' + nevents + '\n'
    fname = f'replay_{i}'
    #print(replay_cmd)
    sbatch_file = open(sbatch_DIR + fname,'w')
    if i == 0 or i == (len(runs_to_replay)-1):
        write_header(sbatch_file,job_name=f'replay_{runs_to_replay[i]}')
    else:
        write_header(sbatch_file,job_name=f'replay_{runs_to_replay[i]}',email=False)
    sbatch_file.write(replay_cmd)
    
    if i == 0:
        sbatch_cmd = f'sbatch current/{fname}\n'
        cmd_file = open(sbatch_DIR + 'commands.txt','w')
        cmd_file.write(sbatch_cmd)
    elif i == (len(runs_to_replay)-1):
        sbatch_cmd = f'sbatch current/{fname}'
        cmd_file = open(sbatch_DIR + 'commands.txt','a')
        cmd_file.write(sbatch_cmd)            
    else:
        sbatch_cmd = f'sbatch current/{fname}\n'
        cmd_file = open(sbatch_DIR + 'commands.txt','a')
        cmd_file.write(sbatch_cmd) 

sbatch_file.close()
cmd_file.close()

with open(sbatch_DIR+'commands.txt','r') as r:
    for l in r:
        #print(l.strip())
        os.system(l.strip())