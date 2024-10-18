import os
import shutil
import subprocess
import math as m
import numpy as np
import myfunctions.functions as mf


class gulp_calc:
    def __init__(
        self,
        hasrun=False,
        rundir='gulp',
        calcnm='general',
        runexe='gulp.x',
        ipname='tmp.gin',
        ouname='tmp.got',
        method='single',
        prism=np.array([10.0, 10.0, 10.0, 90.0, 90.0, 90.0]),
        cortyp='fractional',
        atspec=['H'],
        atspos=np.array([[0.0, 0.0, 0.0]]),
        dpflnm='dump.res',
        maxcyc=1000,
        cutd=3.0,
        cuts=0.8,
        temper=298.15,
        frlbfl='force.lib',
        dmflnm='tmp.res'

        ):

        self.hasrun=hasrun
        self.rundir=rundir
        self.calcnm=calcnm
        self.runexe=runexe
        self.ipname=ipname
        self.ouname=ouname
        self.method=method
        self.prism=prism
        self.cortyp=cortyp
        self.atspec=atspec
        self.atspos=atspos
        self.dpflnm=dpflnm
        self.maxcyc=maxcyc
        self.cutd=cutd
        self.cuts=cuts
        self.temper=temper
        self.frlbfl=frlbfl
        self.dmflnm=self.rundir+'/'+dmflnm

    def print_file(self):
        line = '#comment\n'
#        line = self.print_info()+'\n'
        line = line+self.method+'\n'
        line = line+'name '+self.calcnm+'\n'
        line = line+'cell \n'
        line = line+str(self.prism[0])+' '+str(self.prism[1])+' '+str(self.prism[2])+' '+str(self.prism[3])+' '+str(self.prism[4])+' '+str(self.prism[5])+'\n'
        line = line+self.cortyp+'\n'
        for ii in range(0, len(self.atspec)):
            line = line+self.atspec[ii]+' '+str(self.atspos[ii][0])+' '+str(self.atspos[ii][1])+' '+str(self.atspos[ii][2])+'\n'
        line = line+'maxcyc '+str(self.maxcyc)+'\n'
        line = line+'cutd\n'
        line = line+str(self.cutd)+'\n'
        line = line+'temprature '+str(self.temper)+' K \n'
        line = line+'library nodump '+self.frlbfl+'\n'
        line = line+'dump every '+self.dmflnm+'\n'
        return line

    def run_calc(self):
        os.mkdir(self.rundir)
        fd = open(self.rundir+"/"+self.ipname,'w')
        fd.write(self.print_file())
        fd.close()
        runcmd=self.runexe+' < '+self.rundir+'/'+self.ipname+' > '+self.rundir+'/'+self.ouname
        subprocess.run(runcmd, shell=True)
        self.hasrun=True
        output = self.rundir+'/'+self.ouname
        return

    def get_info(self):
        self.fatpos = []
        self.fatspec = []
        self.fcell = []
        inp = self.rundir+'/'+self.ouname
        print(inp)
        sstring=[]
        sstring.append('Final fractional coordinates of atoms :')
        sstring.append('Final Cartesian lattice vectors (Angstroms) :')
        fd = open(inp, "r").readlines()
        for iline, nline in enumerate(fd):
            if sstring[0] in nline:
                pos_st = iline + 6
            if sstring[1] in nline:
                pos_en = iline - 2
                cel_st = iline + 2
                cel_en = cel_st + 3
        for jj in range(pos_st, pos_en):
            line = fd[jj].split()
            self.fatspec.append(line[1])
            self.fatpos.append(line[3:6])
        self.fatpos = np.array(self.fatpos,dtype=float)
        for jj in range(cel_st,cel_en):
            line = fd[jj].split()
            self.fcell.append(line[0:3])
        self.fcell = np.array(self.fcell,dtype=float)
        self.fprism = cell_convert(self.fcell)
        return
