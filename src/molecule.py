import os
import shutil
import subprocess
import math as m
import random as r
import numpy as np
import myfunctions.functions as mf


class molec:
    def __init__(self,
            name='General',
            #the unit vectors, positions are defined around these vectors so the I matrix just means cartesian.
            cell=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
            #is the system periodic or not
            peri=False,
            #a list of the species present
            spls=[],
            #list of atomic numbers, optional !!I FORGOT HOW I WAS GOING TO IMPLIMENT THESE, PROBABLY ALONG SIDE SPECIES LIST
            atnm=[],
            #list of species masses
            mass=[],
            #index list defining what atom is what species, these integers must match those in spls
            indx=[],
            #atomic positions in terms of the vects
            posi=[],
            #molecular charge
            chrg=0.0,
            #unpaired electrons
            upel=0.0
            ):
        self.name=str(name)
        self.cell=np.array(cell,dtype=float).reshape(3,3)
        self.peri=bool(peri)
        self.spls=np.array(spls).reshape(-1)
        self.atnm=np.array(atnm)
        self.mass=np.array(mass,dtype=float).reshape(-1)
        self.indx=np.array(indx,dtype=int).reshape(-1)
        self.posi=np.array(posi,dtype=float).reshape(-1,3)
        self.chrg=float(chrg)
        self.upel=float(upel)
        return
    #make a copy
    def cpy(self):
        return molec(name=self.name,cell=self.cell,peri=self.peri,spls=self.spls,atnm=self.atnm,mass=self.mass,indx=self.indx,posi=self.posi,chrg=self.chrg,upel=self.upel)
    #shift the molecules
    def shf(self,Sv):
        for ii in range(0,len(self.posi)):
            self.posi[ii]=self.posi[ii]+np.array(Sv)
        return self
    #rotate the molecule
    def rot(self,Rm):
        for ii in range(0,len(self.posi)):
            self.posi[ii]=np.matmul(Rm,self.posi[ii])
        return self
    #replicates a molecule, the vectors can be set OR it will assume the unit vectors should be used
    def repli(self,unit=[1, 1, 1]):
        u=unit
        u=np.array([[u[0],0,0],[0,u[1],0],[0,0,u[2]]],dtype=int) #replications in matrix form

        #make new positions, keep in fractional co-ords, makes things easier
        for ii in range(0,3):
            tmp=self.cpy()
            for jj in range(1,u[ii][ii]):
                vec=np.array(jj*np.identity(3)[ii],dtype=float)
                self.combo(tmp.cpy().shf(vec))
        #make new cell
        u=np.array(u,dtype=float)
        ncell=np.matmul(u,self.cell)
        self.re_cell(ncell)

        return self

    def jiggle(self, inp):
        fac=float(inp)
        for ii in range(0,len(self.posi)):
            for jj in range(0,len(self.posi[ii])):
                off=(r.random()-0.50000)*fac
                self.posi[ii][jj]=self.posi[ii][jj]+off
        return self

    #combine with another molecule
    def combo(self,inpt):
        addt=inpt.cpy()
        self.indx=np.append(self.indx,addt.indx+(len(self.spls)))
        self.spls=np.append(self.spls,addt.spls)
        self.atnm=np.append(self.atnm,addt.atnm)
        self.mass=np.append(self.mass,addt.mass)
        self.posi=np.append(self.posi,addt.posi,axis=0)
        self.chrg=self.chrg+addt.chrg
        self.upel=self.upel+addt.upel
        return self

    #basic print
    def print_bas(self):
        print('structure name: ')
        print(self.name)
        print('basis vectors: ')
        print(self.cell)
        print('prism: ')
        print(mf.cell2prism(self.cell))
        print('periodicity: ')
        print(self.peri)
        print('species: ')
        print(self.spls)
        print('atomic numbers: ')
        print(self.atnm)
        print('masses: ')
        print(self.mass)
        print('charge: ')
        print(self.chrg)
        print('unpaired electrons: ')
        print(self.upel)
        print('multiplicity: ')
        print(self.upel+1)
        print('indexes: ')
        print(self.indx)
        print('positions: ')
        print(self.posi)
        return

    #print structure into xyz file
    def print_xyz(self,fd):
        fd.write(str(len(self.indx))+'\n')
        if self.peri:
            line='Lattice="'
            for ii in self.cell:
                for jj in ii:
                    line=line+str(jj)+' '
            line=line+'" '+self.name+'\n'
            fd.write(line)
        else:
            line=self.name+'\n'
            fd.write(line)

        vec=mf.frac2cart(self.posi,self.cell)
        for ii in range(0,len(self.indx)):
            line=''
            line=line+self.spls[self.indx[ii]]+' '
            for jj in range(0,3):
                line=line+str(vec[ii][jj])+' '
            line=line+'\n'
            fd.write(line)
        return

    #print the structure into a gen file
    def print_gen(self,fd):
        if self.peri:
            fd.write(str(len(self.indx))+' F \n')

            line=''
            for ii in self.spls:
                line=line+ii+' '
            line=line+'\n'
            fd.write(line)

            for ii in range(0,len(self.indx)):
                line=''
                line=line+str(ii+1)+' '
                line=line+str(self.indx[ii]+1)+' '
                for jj in range(0,3):
                    line=line+str(self.posi[ii][jj])+' '
                line=line+'\n'
                fd.write(line)

            line='1.0 1.0 1.0 \n'
            for jj in range(0,3):
                for kk in range(0,3):
                    line=line+str(self.cell[jj][kk])+' '
                line=line+'\n'
            fd.write(line)

        else:
            fd.write(str(len(self.indx))+' C \n')

            line=''
            for ii in self.spls:
                line=line+ii+' '
            line=line+'\n'
            fd.write(line)
            for ii in range(0,len(self.indx)):
                vec=np.matmul(self.cell,self.posi[ii])
                line=''
                line=line+str(ii+1)+' '
                line=line+str(self.indx[ii]+1)+' '
                for jj in range(0,3):
                    line=line+str(vec[jj])+' '
                line=line+'\n'
                fd.write(line)
        return
    #only keep the atoms specified in the provided vecor
    def keep(self, vec):
        vec=np.array(vec,dtype=int)
        newpos=[]
        newind=[]
        for ii in vec:
            newpos.append(self.posi[ii])
            newind.append(self.indx[ii])
        #will also need to keep masses and charges when i impliment them properly
        self.indx=np.array(newind,dtype=int).reshape(-1)
        self.posi=np.array(newpos,dtype=float).reshape(-1,3)
        return self

    def remove(self, vec):
        k=[]
        for ii in range(0,len(self.posi)):
            if ii not in vec:
                k.append(ii)
        k=np.array(k,dtype=int)
        return self.keep(k)

    #shrink the species and index list so they a minimal
    def shrink(self):
        #first, get shrunk species list
        nlst=[]
        nlst.append(self.spls[0])
        for ii in self.spls:
            add=True
            for jj in nlst:
                if ii == jj:
                    add=False
            if add:
                nlst.append(ii)


        nind=[]
        #now get the index
        for ii in range(0,len(self.indx)):
            for jj in range(0,len(nlst)):
                if self.spls[self.indx[ii]] == nlst[jj]:
                    nind.append(jj)

        self.spls=np.array(nlst)
        self.indx=np.array(nind,dtype=int)
        return self

    #function for giving the structure a new cell wil maintaining atomic cartesian positions
    def re_cell(self,ncell):
        #convert from fractional to cartesian
        cart=mf.frac2cart(self.posi,self.cell)
        #swap cells
        self.cell=np.array(ncell,dtype=float).reshape(3,3)
        #convert positions back
        self.posi=mf.cart2frac(cart,self.cell)

        return self

    #function for changing our cell allignment and making slabs
    #only works on positive indicies, needs work on formalising the math and cleaning up, but it "works" right now
    #absolute garbage but abbi was sick.
    def realign(self,miller=[1, 1, 1]):
        u=np.array(miller,dtype=int)
        #very similar to our replicator but lets bring the tmp definition outside the first loop
        tmp=self.cpy()
        #first, get length of final z vector
        lzvec=np.linalg.norm(miller@self.cell)
        #now i want self to be blank, just the unit cell
        self=molec(name=self.name,cell=self.cell,peri=self.peri)
        vec=np.array([0.0, 0.0, 0.0],dtype=float)
        for ii in range(0,3):
            for jj in range(0,u[ii]):
                self.combo(tmp.cpy().shf(vec))
                vec=vec+np.array(np.identity(3)[ii],dtype=float)
        #align main vector with the Z axis
        zvec=np.array([0.0, 0.0, 1.0])
        avec=np.cross(zvec,vec)
        thet=np.rad2deg(m.asin(np.linalg.norm(avec)/(np.linalg.norm(vec)*np.linalg.norm(zvec))))
        Rm=mf.Rg(avec,thet)
        self.cell=self.cell@Rm
        #align b vector with y direction
        yvec=np.array([0.0, 1.0, 0.0])
        bvec=np.array(self.cell[1])
        bvec[2]=0
        avec=np.cross(yvec,bvec)
        thet=np.rad2deg(m.asin(np.linalg.norm(avec)/(np.linalg.norm(bvec)*np.linalg.norm(yvec))))
        Rm=mf.Rg(zvec,thet)
        self.cell=self.cell@Rm
        #make new cell
        mmat=np.array([[0.0, u[2], -1*u[1]],[-1*u[2], 0.0, u[0]],[u[1], -1*u[0], 0.0]])
        ncell=np.matmul(mmat,self.cell)
        #make candidate cells
        ccell=[]
        for ii in [[0,1], [0,2], [1,2]]:
            ccell.append(np.array([ncell[ii[0]], ncell[ii[1]], lzvec*zvec]))
        
        for ii in ccell:
            if np.linalg.det(ii) < 0.0:
                ii[0]=-1.0*ii[0]
            if np.linalg.det(ii) > 1.0e-10:
                ncell=ii
        self.re_cell(ncell)
        return self

#function reads a pdb file and returns a list or structures, will need modifying but lets start with something simple for now
def pdb_read(fd):
    #keywords we are looking for
    cstr='CRYST1'
    astr='ATOM'
    estr='END'

    lines = fd.readlines()

    #the structures
    strucs=[]

    pvec=[] #unit cell in prism form
    atlst=[] #atomic positions
    slst=[] #species list
    for upto in range(0,len(lines)):
        #look for unit cell definition
        if lines[upto].find(cstr)!= -1:
            pvec=np.array(lines[upto].split()[1:7]).astype(float)
        #find all the atoms
        if lines[upto].find(astr)!= -1:
            atlst.append(np.array(lines[upto].split()[5:8]).astype(float))
            slst.append(lines[upto].split()[10])
        #if you reach an end, make a molecule and reset
        if lines[upto].find(estr)!= -1:
            strucs.append(molec(peri=True,cell=mf.prism2cell(pvec),spls=slst,indx=range(0,len(slst)),posi=mf.cart2frac(atlst,mf.prism2cell(pvec))).shrink())
            pvec=[]
            atlst=[]
            slst=[]

    return strucs

#function reads a gin file and returns the molecular structure, dont really want to use it right now
#def gin_read(fd):


#read in an xyz file into a list of structures
def read_xyz(fd):
    #keywords we are looking for

    lines = fd.readlines()

    #the structures
    strucs=[]

    cell=[] #unit cell
    posi=[] #atomic positions
    spls=[] #species list
    upto=0
    while upto < len(lines):
        #get the number of atoms in this structure
        atnum=int(lines[upto].split()[0])
        #skip comment
        upto=upto+2
        #find all the atoms
        for ii in range(upto,upto+atnum):
            posi.append(np.array(lines[ii].split()[1:4]).astype(float))
            spls.append(lines[ii].split()[0])
        #if you reach an end, make a molecule and reset
        inpt=molec(peri=False,spls=spls,indx=range(0,len(spls)),posi=posi).shrink()
        strucs.append(inpt)
        cell=[]
        posi=[]
        spls=[]
        upto=upto+atnum

    return strucs

