import os
import shutil
import subprocess
import math as m
import numpy as np
import myfunctions.functions as mf
import myfunctions.molecule as ml

class calc:
    def __init__(
        self,
        moduse=os.environ.get('LOADEDMODULES').split(':'),
        hasrun=False,
        rundir='/.',
        runexe=os.environ.get('CP2K_BASE')+'/bin/cp2k.psmp',
        ipname='cp2k.in',
        ouname='cp2k.out',
        
        #input structure
        inmole=ml.molec(name='general'),

        #global block options
        #because there is a bunch of redundancy between methods i think its a good idea to just collect variables and then use other means to organise the output file

        runtyp='ENERGY_FORCE',
        projec='general',
        prntlv='LOW',
        global_sec=True,
        
        #mpi parameter
        mpiexe=os.environ.get('OMPI_BASE')+'/bin/mpirun',
        ucores=32,
        
        #files and parameters for dealing with submitting caclualtions
        runtim='24:00:00',
        usemem='100G',
        finfln='cp2k.fin',
        subfln='cp2k.sh'
        ):
        self.hasrun=hasrun
        self.moduse=moduse
        self.rundir=os.environ.get('PWD')+'/'+rundir
        self.runexe=runexe
        self.ipname=self.rundir+'/'+ipname
        self.ouname=self.rundir+'/'+ouname
            
        self.inmole=inmole
        #self.kponts=np.array(kponts,dtype=float).reshape(-1,4)
        self.kponts=kponts

        self.control_sec=control_sec
        self.calcul=calcul
        self.restar=restar
        self.prefix=prefix
        self.pseudo=pseudo
        self.outdir=self.rundir+'/'+outdir
        self.verbos=verbos

        self.system_sec=system_sec
        self.ibravs=ibravs
        self.noncol=noncol
        self.lspino=lspino
        self.nosyme=nosyme
        self.ecutwf=ecutwf
        self.ecutrh=ecutrh
        self.occupa=occupa
        self.smeari=smeari
        self.degaus=degaus

        self.electron_sec=electron_sec
        self.elemax=elemax
        self.conthr=conthr
        self.mixbet=mixbet

        if self.calcul in ['relax', 'md', 'vc-relax', 'vc-md']:
            self.ion_sec=True
        else:
            self.ion_sec=ion_sec
        self.iondyn=iondyn

        if self.calcul in ['vc-relax', 'vc-md'] :
            self.cell_sec=True
        else:
            self.cell_sec=cell_sec
        self.celdyn=celdyn
        self.spepse=spepse
            
        self.mpiexe=mpiexe
        self.ucores=ucores
        self.runtim=runtim
        self.usemem=usemem

        self.finfln=self.rundir+'/'+finfln
        self.subfln=self.rundir+'/'+subfln

        return

    def print_bas(self):
        return self

    def print_inpt(self,fd):
        #make the input file, for quantum_E its all one file
        if self.control_sec:
            fd.write('&control \n')
            fd.write(' calculation = \''+self.calcul+'\' \n')
            fd.write(' restart_mode = \''+self.restar+'\' \n')
            fd.write(' prefix = \''+self.prefix+'\' \n')
            fd.write(' pseudo_dir = \''+self.pseudo+'\' \n')
            fd.write(' outdir = \''+self.outdir+'\' \n')
            fd.write(' verbosity = \''+self.verbos+'\' \n')
            fd.write('/\n')
        if self.system_sec:
            fd.write('&system \n')
            fd.write(' ibrav = '+str(self.ibravs)+' \n')
            fd.write(' noncolin = '+str(self.noncol)+' \n')
            fd.write(' lspinorb = '+str(self.lspino)+' \n')
            fd.write(' nosym = '+str(self.nosyme)+' \n')
            fd.write(' nat = '+str(len(self.inmole.indx))+' \n')
            fd.write(' ntyp = '+str(len(self.inmole.spls))+' \n')
            fd.write(' ecutwfc = '+str(self.ecutwf)+' \n')
            fd.write(' ecutrho = '+str(self.ecutrh)+' \n')
            fd.write(' occupations = \''+self.occupa+'\' \n')
            fd.write(' smearing = \''+self.smeari+'\' \n')
            fd.write(' degauss = '+str(self.degaus)+' \n')
            fd.write('/\n')
        if self.electron_sec:
            fd.write('&electrons \n')
            fd.write(' electron_maxstep = '+str(self.elemax)+' \n')
            fd.write(' conv_thr =  '+str(self.conthr)+' \n')
            fd.write(' mixing_beta = '+str(self.mixbet)+' \n')
            fd.write('/\n')
        if self.ion_sec:
            fd.write('&ions \n')
            fd.write(' ion_dynamics = \''+self.iondyn+'\' \n')
            fd.write('/\n')
        if self.cell_sec: 
            fd.write('&cell \n')
            fd.write(' cell_dynamics = \''+self.celdyn+'\' \n')
            fd.write('/\n')
        fd.write('ATOMIC_SPECIES \n')
        for ii in self.inmole.spls:
            fd.write(' '+ii+' '+str(self.spepse[ii][0])+' '+self.spepse[ii][1]+' \n')
        fd.write('ATOMIC_POSITIONS {alat} \n')
        for ii in range(0,len(self.inmole.indx)):
            line=' '+self.inmole.spls[self.inmole.indx[ii]]+' '
            for jj in self.inmole.posi[ii]:
                line=line+str(jj)+' '
            line=line+' \n'
            fd.write(line)
        fd.write('CELL_PARAMETERS {angstrom} \n')
        for ii in self.inmole.cell:
            line=' '
            for jj in ii:
                line=line+str(jj)+' '
            line=line+'\n'
            fd.write(line)
        fd.write(self.kponts) 
        #eventually ill use my own set of kpoint generators but right now its just going to be text
        #fd.write('K_POINTS {tpiba} \n')
        #fd.write(' '+str(len(self.kponts))+' \n')
        #for ii in self.kponts:
        #    line=' '
        #    for jj in ii:
        #        line=line+str(jj)+' '
        #    line=line+'\n'
        #    fd.write(line)
        return self

    def print_subfl(self,fd):
        line='#!/bin/bash\n'
        line=line+'#PBS -P f97\n'
        line=line+'#PBS -N '+self.inmole.name+'_'+self.prefix+'\n'
        line=line+'#PBS -l walltime={0}\n'.format(self.runtim)
        line=line+'#PBS -l ncpus='+str(self.ucores)+'\n'
        line=line+'#PBS -l mem='+self.usemem+'\n'
        line=line+'#PBS -l wd\n'
        line=line+'#PBS -l storage=gdata/f97+scratch/f97\n'
        line=line+'module load '
        for ii in self.moduse:
            line=line+ii+' '
        line=line+'\n'
        line=line+'module list\n'
        line=line+self.mpiexe+' -c ${PBS_NCPUS} '+self.runexe+' < '+self.ipname+' > '+self.ouname+'\n'
        line=line+'touch '+self.finfln+'\n'
        fd.write(line)
        return self


    def run_calc(self):
        return self

    def set_que_calc(self):
        os.mkdir(self.rundir)
        fd = open(self.ipname,'w')
        self.print_inpt(fd)
        fd.close()
        fd = open(self.subfln,'w')
        self.print_subfl(fd)
        fd.close()
        return self


    def que_calc(self):
        maindir=os.getcwd()
        self.set_que_calc()
        runcmd='cd '+self.rundir+' ; qsub '+self.subfln+' ; cd '+maindir
        subprocess.run(runcmd, shell=True)
        self.hasrun=True
        return

    def get_info(self):
        return

