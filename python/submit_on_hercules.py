#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
import subprocess

currentDir = os.getcwd();
CMSSWDir = currentDir+"../";
ReducedTreeDir = "";

name = ["WJets","TTbar","RSGraviton1000"];

for i in range(len(name)):
    fn = "Job/Job_"+name[i];
    outScript = open(fn+".sh","w");
    command = "python python/produceWWNtuples.py -i /gwteray/users/brianza/WWNtupleRun2/"+name[i]+".root -o WWTree_"+name[i];
    outScript.write('#!/bin/bash');
    outScript.write("\n"+'cd '+CMSSWDir);
    outScript.write("\n"+'eval `scram runtime -sh`');
    outScript.write("\n"+'cd '+currentDir);
    outScript.write("\n"+"unbuffer "+command+" > "+fn+"_output.txt");
    outScript.close();
    os.system("chmod 777 "+currentDir+"/"+fn+".sh");
    os.system("qsub -V -d "+currentDir+" -q shortcms "+currentDir+"/"+fn+".sh");