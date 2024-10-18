import os
import shutil
import subprocess
import math as m
import numpy as np
import myfunctions.functions as mf
import myfunctions.molecule as ml

class sies_calc:
    def __init__(
            self,
            hasrun=False,
            rundir='sies',
            calcnm='general',
            callbl='general',
            runexe='siesta',
            ipname='sies.in',
            ouname='sies.out',
            ername='sies.err'
            #input structure
            inmole=ml.molec()
            #mpi parameter
            mpiexe='mpirun',
            ucores=1,
            #files and parameters for dealing with submitting caclualtions
            runtim=[24, 0, 0],
            usemem=100,
            finfln='fin.sies',
            subfln='sies.sh'
            ):
        self.hasrun=hasrun
        self.rundir=rundir
        self.calcnm=calcnm
        self.runexe=runexe
        self.ipname=ipname
        self.ouname=ouname
        self.inmole=inmole
        self.mpiexe=mpiexe
        self.ucores=ucores
        self.usemem=usemem
        self.finfln=finfln
        self.subfln=subfln
        return

    def print_inpt(self):
        maindir=os.getcwd()
        os.mkdir(self.rundir)
        
        
        fd = open(self.rundir+'/'+self.ipname,'w')

        fd.close()
        return

    def print_sbfl(self):
        line='#!/bin/bash\n'
        line=line+'#PBS -P f97\n'
        line=line+'#PBS -N '+self.calcnm+'\n'
        line=line+'#PBS -l walltime=24:00:00\n'
        line=line+'#PBS -l ncpus='+str(self.cores)+'\n'
        line=line+'#PBS -l mem=100GB\n'
        line=line+'#PBS -l wd\n'
        line=line+'#PBS -l storage=gdata/f97+scratch/f97\n'
        line=line+'#PBS -o '+self.ouname+'\n'
        line=line+'#PBS -e dftb.err\n'
        line=line+'module load openmpi/4.1.2 dftbplus/22.1\n'
        line=line+'module list\n'
        line=line+self.mpiexe+' -c ${PBS_NCPUS} '+self.runexe+'\n'
        line=line+'touch '+self.finflnm+'\n'
        return line


    def run_calc(self):
        maindir=os.getcwd()
        os.mkdir(self.rundir)
        fd = open(self.rundir+'/'+self.ipname,'w')
        fd.write(self.print_file())
        fd.close()
        runcmd='cd '+self.rundir+' ; '+self.runexe+' > '+self.ouname+' ; touch '+self.finflnm+'; cd '+maindir
        subprocess.run(runcmd, shell=True)
        self.hasrun=True
        return

    def que_calc(self):
        maindir=os.getcwd()
        os.mkdir(self.rundir)
        fd = open(self.rundir+'/'+self.ipname,'w')
        fd.write(self.print_file())
        fd.close()
        fd = open(self.rundir+'/'+self.subflnm,'w')
        fd.write(self.print_subfl())
        fd.close()
        runcmd='cd '+self.rundir+' ; qsub '+self.subflnm+' ; cd '+maindir
        subprocess.run(runcmd, shell=True)
        self.hasrun=True
        return

    def get_info(self):
        return
