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
            hasrun=False,
            rundir='dftb',
            calcnm='general',
            callbl='general',
            runexe='dftb+',
            ipname='dftb_in.hsd',
            ouname='dftb.out',
            ername='dftb.err',
            outprf='geo_end.gen',
            #input structure
            inmole=ml.molec(),
            #driver options
            driver='CongugateGradient',
            #calculations parameters
            kponts=np.array([[1, 0, 0],[0, 1, 0],[0, 0, 1]]),
            #mpi parameter
            mpiexe='mpirun',
            ucores=1,
            #files and parameters for dealing with submitting caclualtions
            runtim=[24, 0, 0],
            usemem=100,
            finfln='dftb.fin',
            subfln='dftb.sh'
            ):
        self.hasrun=hasrun
        self.rundir=rundir
        self.calcnm=calcnm
        self.callbl=callbl
        self.runexe=runexe
        self.ipname=ipname
        self.ouname=ouname
        self.ername=ername
        self.outprf=outprf
        self.inmole=inmole
        self.driver=driver
        self.kponts=kponts
        self.mpiexe=mpiexe
        self.ucores=ucores
        self.runtim=runtim
        self.usemem=usemem
        self.finfln=finfln
        self.subfln=subfln
        return

    def print_bas(self):
        return

    def print_inpt(self,fd):

        #write geometry input
        self.print_geomet(fd)
        #write geometry driver
        self.print_driver(fd)
        #write hamiltonian parameters
        #self.print_hamilt(fd)
        #

        return

    #print geometry part of input file
    def print_geomet(self,fd):
        fd.write("Geometry = GenFormat {\n")
        self.inmole.print_gen(fd)
        fd.write("}\n")
        return

    def print_driver(self,fd):
        fd.write("Driver = "+self.driver+" { \n")
        if self.driver == "CongugateGradient":
            fd.write("  AppendGeometries = Yes\n")
            fd.write("  LatticeOpt = Yes\n")
            fd.write("  MaxForceComponent = 1.0d-04\n")
            fd.write("  MaxSteps = 100000\n")
            fd.write("  MovedAtoms = 1:-1\n")
            fd.write("  MaxAtomStep = 0.200000000000000\n")
            fd.write("  StepSize = 100.000000000000\n")
            fd.write("  OutputPrefix = 'latopt'\n")
            fd.write("  Constraints = {}\n")
            fd.write("  Pressure = 0.000000000000000E+000\n")
            fd.write("  FixAngles = No\n")
            fd.write("  Isotropic = No\n")
            fd.write("  MaxLatticeStep = 0.200000000000000\n")
        elif self.driver == "VelocityVerlet":
            fd.write("  TimeStep [Femtosecond] = 1.0\n")
            fd.write("  MovedAtoms = 1:-1\n")
            fd.write("  MDRestartFrequency = 10\n")
            fd.write("  Thermostat = Berendsen {\n")
            fd.write("    TimeScale [fs] = 100\n")
            fd.write("    AdaptFillingTemp = Yes\n")
            fd.write("    Temperature [Kelvin] = TemperatureProfile {\n")
            fd.write("constant     1000000  300.0\n")
            fd.write("    }\n")
            fd.write("  }\n")
            fd.write("  Velocities = {}\n")
            fd.write("  KeepStationary = Yes\n")
            fd.write("  OutputPrefix = 'geo_end'\n")
            fd.write("  Plumed = No\n")
            fd.write("  Barostat = {\n")
            fd.write("    Pressure [Pa] = 100000.0\n")
            fd.write("    Timescale [fs] = 100\n")
            fd.write("    Isotropic = Yes\n")
            fd.write("  }\n")
            fd.write("  Masses = {}\n")
        return

    def print_subfl(self):
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
