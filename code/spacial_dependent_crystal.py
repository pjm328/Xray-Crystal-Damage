from math import *
import os
import multiprocessing as mp

mlfsom_path='/home/peter/Desktop/Summer_Research_2019/mlfsom/'
code_path='/home/peter/Desktop/Summer_Research_2019/code/'
def gaussian(peak):
    return lambda x,y: peak*e**(-x**2-y**2)

def get_experiment_list(N,start_mos,k_mos,k_cell,k_bfactor,beam_func,frames,frame_rate,radial_symmetry=True,
    angle_slices=1,osc=0.01):
    '''Returns a queue of parameters of experiments, weights on the diffraction pattern,
    and frame number, for mlfsom to calculate the spacial and dose
    dependent diffraction patterns for a crystal model given by:
    N by N grid of cells that are centered on the beam
    k_mos - increase in mosaicity given by mos_0 + k_mos*dose
    k_cell - increase in unit cell dimensions given by cell_dimenion_0 + k_cell*dose
    k_bfactor - increase in b-factors given by B_0 + k_bfactor*dose
    beam_func - defines spacial distribution of power of beam
    radial_symmetry - if True, assumes beam_func is radially symmetric, if False, does not
    frames - number of frames
    frame_rate - frame rate of experiment
    Returns tuple of list of parameters for experiments to be calculated by mlfsom given by (ID,mos,bfactor_inc,cell_inc,frames,osc)
    and list of weights to each set of experiment on each frame given by (frame,weight,exp_ID)
    '''
    #create list of unique cells and their weights in the crystal
    cell_list=[] #form is (x,y,weight)
    if radial_symmetry:
        #less cells are needed due to symmetry
        if N%2==0:
            for i in range(0,N//2):
                for j in range(0,i+1):
                    if i==j:
                        weight=4
                    else:
                        weight=8
                    cell_list.append((i+.5,j+.5,weight))
        else:
            for i in range(0,N//2+1):
                for j in range(0,i+1):
                    if (i,j)==(0,0):
                        weight=1
                    elif i==j or j==0:
                        weight=4
                    else:
                        weight=8
                    cell_list.append((i,j,weight))
    else:
        #all cells are needed
        if N%2==0:
            for i in range(-N//2,N//2):
                for j in range(-N//2,N//2):
                    cell_list.append((i+.5,j+.5,1))
        else:
            for i in range(-(N//2),N//2+1):
                for j in range(-(N//2),N//2+1):
                    cell_list.append((i,j,1))
    #calculate dose in each cell for each frame and parameters for each cell
    experiment_list=[] #form is (ID,mos,bfactor_inc,cell_inc,frames,osc)
    frame_weights=[] #form is (frame,weight,experiment_ID)
    for frame in range(0,frames):
        exposure=frame/frame_rate
        for cell in cell_list:
            ID=len(experiment_list)
            frame_weights.append((frame+1,cell[2],ID))
            dose=beam_func(cell[0],cell[1])*exposure
            experiment_list.append((ID,start_mos+k_mos*dose,k_bfactor*dose,k_cell*dose,angle_slices,osc))
    #remove repeated experiments and generate list of weights
    i=0
    while i<len(experiment_list):
            experiment=experiment_list[i]
            j=i+1
            while j<len(experiment_list):
                other=experiment_list[j]
                if experiment[1:6]==other[1:6]:
                    experiment_list.pop(j)
                    for index in range(len(frame_weights)):
                        params=frame_weights[index]
                        if params[2]==other[0]:
                            frame_weights[index]=(params[0],params[1],experiment[0])
                else:
                    j+=1
            i+=1
    return (frame_weights,experiment_list)

def get_homogenous_experiment_list(start_mos,k_mos,k_cell,k_bfactor,frames,frame_rate,angle_slices=1,osc=0.01):
    experiment_list=[] #form is (ID,mos,bfactor_inc,cell_inc,frames,osc)
    for frame in range(0,frames):
        exposure=frame*frame_rate
        for cell in cell_list:
            ID=len(experiment_list)
            dose=1*exposure
            experiment_list.append((ID,start_mos+k_mos*dose,k_bfactor*dose,k_cell*dose,angle_slices,osc))
    return experiment_list

def run_experiments(experiment_list,prefix,threads):
    #copy pdb file to temp file
    prev_cwd=os.getcwd()
    os.chdir(mlfsom_path)
    os.system('cp 1H87.pdb temp.pdb')
    #clear previous temp files
    os.system('rm /tmp/peter/*')
    #set up threading and run experiments
    p=mp.Pool(threads)
    experiments_and_prefix=list(zip(experiment_list,[prefix]*len(experiment_list)))
    for i in range(0,len(experiments_and_prefix),threads):
        #generate input files
        for exp in experiments_and_prefix[i:i+threads]:
            ID=exp[0][0]
            bfactor_inc=exp[0][2]
            cell_inc=exp[0][3]
            print("generating inputs for id="+str(ID))
            os.system('./change_pdb_param.com temp.pdb input'+str(ID)+'.pdb add_cell='+str(cell_inc)+
            ' add_bfactor='+str(bfactor_inc))
            os.system('./ano_sfall.com energy=12660 input'+str(ID)+'.pdb 1.2A solvent_B=35')
            os.system('mv ideal_ano.mtz input'+str(ID)+'.mtz')
        p.map(run_experiment,experiments_and_prefix[i:i+threads])
        #clear temp files and move results
        print("clearing temp files and moving results")
        os.system('mv /tmp/peter/mlfsom*.XYI ~/Desktop/Summer_Research_2019/mlfsom/data/tempdata')
        os.system('mv /tmp/peter/mlfsom*predin.txt ~/Desktop/Summer_Research_2019/mlfsom/data/tempdata')
        os.system('rm ~/Desktop/Summer_Research_2019/mlfsom/fit2d_*')
        os.system('rm /tmp/peter/*')
    os.chdir(prev_cwd)

def run_experiment(experiment_and_prefix):
    experiment=experiment_and_prefix[0]
    prefix=experiment_and_prefix[1]
    ID=experiment[0]
    mos=experiment[1]
    bfactor_inc=experiment[2]
    cell_inc=experiment[3]
    frames=experiment[4]
    osc=experiment[5]
    print('running id='+str(ID)+' mos='+str(mos)+' bfactor_inc='+str(bfactor_inc)+' cell_inc='+str(cell_inc))
    os.system('./mlfsom.com ~/Desktop/Summer_Research_2019/mlfsom/data/'+prefix+str(ID)+
    '_001.img input'+str(ID)+'.mtz mosaic='+str(mos)+' frames='+str(frames)+' osc='+str(osc))
    print('----id: '+str(ID)+' is done----')

def spacial_dependent_crystal(N,start_mos,k_mos,k_cell,k_bfactor,beam_func,frames,frame_rate,prefix,threads=3,angle_slices=1,osc=.01,radial_symmetry=True):
    frame_weights,experiment_list=get_experiment_list(N,start_mos,k_mos,k_cell,k_bfactor,beam_func,frames,frame_rate,radial_symmetry=radial_symmetry,
    angle_slices=angle_slices,osc=osc)
    write_description(N,start_mos,k_mos,k_cell,k_bfactor,beam_func,frames,frame_rate,prefix,angle_slices,osc,frame_weights)
    write_exp_queue(experiment_list,prefix)
    #run_experiments(experiment_list,prefix,threads)
    return experiment_list

def homogenous_crystal(start_mos,k_mos,k_cell,k_bfactor,frames,frame_rate,prefix,threads=3,angle_slices=1,osc=.01):
    experiment_list=get_homogenous_experiment_list(start_mos,k_mos,k_cell,k_bfactor,frames,frame_rate,angle_slices=angle_slices,osc=osc)
    run_experiments(experiment_list,prefix,threads)
    return experiment_list

def write_exp_queue(experiment_list,prefix):
    f=open("exp-"+prefix+"-queue.txt","w")
    for exp in experiment_list:
        f.write(str(exp)+"\n")
    f.close()

def run_from_exp_queue(queue_file,N,threads,prefix):
    exp_list=[]
    with open(queue_file,"r") as f, open("tmp"+queue_file,"w") as out:
        for i in range(N):
            exp_list.append(eval(next(f)))
        for line in f:
            out.write(line)
    os.remove(queue_file)
    os.rename("tmp"+queue_file,queue_file)
    run_experiments(exp_list,prefix,threads)

def write_description(N,start_mos,k_mos,k_cell,k_bfactor,beam_func,frames,frame_rate,prefix,angle_slices,osc,frame_weights):
    f=open("exp-"+prefix+"-desc.txt","w")
    f.write("experiment - "+prefix+"\n")
    f.write("N="+str(N)+"\n")
    f.write("Start Mos="+str(start_mos)+"\n")
    f.write("k_mos="+str(k_mos)+"\n")
    f.write("k_cell="+str(k_cell)+"\n")
    f.write("k_bfactor="+str(k_bfactor)+"\n")
    f.write("beam_func="+str(beam_func)+"\n")
    f.write("frames="+str(frames)+"\n")
    f.write("frame_rate="+str(frame_rate)+"\n")
    f.write("angle_slices="+str(angle_slices)+"\n")
    f.write("osc="+str(osc)+"\n")
    f.write("exp list weights:"+"\n")
    f.write("(frame,weight,experiment_ID)"+"\n")
    for i in frame_weights:
        f.write(str(i)+"\n")
    f.close()

spacial_dependent_crystal(5,0.05,0.07,0.0153,2.228,gaussian(10/14),14,1,"sp_seq1_")
