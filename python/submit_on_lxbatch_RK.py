#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
import subprocess

currentDir = os.getcwd();
CMSSWDir = currentDir+"/../";
ReducedTreeDir = "";

name = ["RSGraviton600", "RSGraviton800", "RSGraviton1000", "RSGraviton1200", "RSGraviton1400", "RSGraviton1600", "RSGraviton1800", "RSGraviton2000",
        "RSGraviton2500", "RSGraviton3000", "RSGraviton3500", "RSGraviton4000", "RSGraviton4500",
        "WJets","WJets100", "WJets200", "WJets400", "WJets600", "TTbar_amcatnlo", "sch", "tch", "tch_bar", "tWch", "tWch_bar", "WW", "ZZ", "WZ",
        "WW_excl", "WZ_excl", "ZZ_excl", "TTbar_powheg", "TTbar_madgraph",
        "WJets_50ns", "WW_50ns", "WZ_50ns", "ZZ_50ns", "TTbar_amcatnlo_50ns", "TTbar_powheg_50ns", "TTbar_madgraph_50ns",
        "tch_50ns", "tch_bar_50ns", "tWch_bar_50ns", "tWch_50ns",
        "BulkGraviton800", "BulkGraviton1000", "BulkGraviton1200", "BulkGraviton1400", "BulkGraviton1600", "BulkGraviton1800", 
        "BulkGraviton2000", "BulkGraviton2500", "BulkGraviton3000", "BulkGraviton4000", "BulkGraviton4500",
        "WprimeToWZ600","WprimeToWZ800","WprimeToWZ1000","WprimeToWZ1200","WprimeToWZ1400","WprimeToWZ1600","WprimeToWZ1800",
        "WprimeToWZ2000","WprimeToWZ2500","WprimeToWZ3000","WprimeToWZ3500","WprimeToWZ4000","WprimeToWZ4500"];

category = ["mu","el"];
xSecWeight = ["4.76735", "1.16691", "0.377865", "0.144482", "0.0616708", "0.0288651", "0.0141334", "0.00751431",
              "0.00167726", "0.000443483", "0.000133915", "0.0000424117", "0.0000130705",
              "61526.7", "1292.", "385.9", "47.9", "19.9", "831.76", "10.11", "43.8", "26.07", "35.6", "35.6", "118.7", "15.4", "66.1",
              "43.53", "10.96", "3.38", "831.76", "831.76",
              "61526.7", "118.7", "15.4", "66.1", "831.76", "831.76", "831.76",
              "43.8", "26.07", "35.6", "35.6",
              "0.001332687", "0.000359194", "0.000119842", "0.000045798", "0.", "0.", 
              "0.000004197", "0.000000786", "0.000000172", "0.", "0.",
              "0.","0.","0.","0.","0.","0.","0.",
              "0.","0.","0.","0.","0.","0."];

N = ["32354.", "31906.", "32448.", "32252.", "32275.", "31971.", "32021.", "31295.",
     "32032.", "31374.", "32194.", "32207.", "31551.",
     "24151270.","10142187.", "5231856.", "1901705.", "1036108.", "42730273.", "984400.", "2966200.", "1695400.", "995600.", "1000000.", "994416.", "996168.", "991232.",
     "1969600.", "24711046.", "18898680.", "19899500.", "11339232.",
     "24089991.", "989608.", "998848.", "996920.", "4994250.", "19665194.", "4992231.",
     "1273800.", "681900.", "1000000.", "998400.",
     "50000.", "50000.", "50000.", "50000.", "49200.", "48400.", 
     "50000.", "48400.", "49800.", "50000.", "50000.",
              "0.","0.","0.","0.","0.","0.","0.",
              "0.","0.","0.","0.","0.","0."];

mass = ["600", "800", "1000", "1200", "1400", "1600", "1800", "2000",
        "2500", "3000", "3500", "4000", "4500",
        "0","0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0",
        "0", "0", "0", "0", "0",
        "0", "0", "0", "0", "0", "0", "0",
        "0", "0", "0", "0",
        "800", "1000", "1200", "1400", "1600", "1800",
        "2000", "2500", "3000", "4000", "4500",
        "600", "800", "1000", "1200", "1400", "1600", "1800",
        "2000", "2500", "3000", "3500", "4000", "4500"];

nameData = ["data_mu_prompt","data_el_prompt","data_mu_17jul","data_el_17jul"];

#MC
for a in range(len(category)):
    for i in range(len(name)):
        fn = "Job/Job_"+name[i]+"_"+category[a];
        outScript = open(fn+".csh","w");
        command = "python python/produceWWNtuples.py -n "+name[i]+" -o WWTree_"+name[i]+" -l "+category[a]+" -w "+xSecWeight[i]+" -no "+N[i]+" -mass "+mass[i]+" --ismc True -trig 1";
        print command;
        outScript.write('#!/bin/bash');
        outScript.write("\n"+'cd '+CMSSWDir);
        outScript.write("\n"+'eval `scram runtime -sh`');
        outScript.write("\n"+'cd '+currentDir);
        outScript.write("\n"+command);
        outScript.close();
        os.system("chmod 777 "+currentDir+"/"+fn+".csh");
        #command2 = "bsub -q 1nd -cwd "+currentDir+" "+currentDir+"/"+fn+".sh";
        #os.system(command2);
        #print command2

for a in range(len(category)):
    for i in range(len(name)):
		fn = "Job/Condor_"+name[i]+"_"+category[a];
		outScript = open(fn,"w");
		outScript.write('universe = vanilla');
		outScript.write('\nExecutable = Job/Job_'+name[i]+'_'+category[a]+'.csh');
		outScript.write('\nRequirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000');
		outScript.write('\nShould_Transfer_Files = YES');
		outScript.write('\nWhenToTransferOutput = ON_EXIT');
		outScript.write('\nOutput		= condor_logs/'+name[i]+'_$(Cluster)_$(Process).stdout');
		outScript.write('\nError		= condor_logs/'+name[i]+'_$(Cluster)_$(Process).stderr');
		outScript.write('\nLog		= condor_logs/'+name[i]+'_$(Cluster)_$(Process).log');
		outScript.write('\nQueue 1');
		outScript.close();
		command2 = "condor_submit "+currentDir+"/"+fn;
		os.system(command2);
		print command2


#data
for a in range(len(category)):
    for i in range(len(nameData)):
        fn = "Job/Job_"+nameData[i]+"_"+category[a];
        outScript = open(fn+".csh","w");
        command = "python python/produceWWNtuples.py -n "+nameData[i]+" -o WWTree_"+nameData[i]+" -l "+category[a]+" -w 1. -no 1. -mass 0 --ismc False -trig 1";
        print command;
        outScript.write('#!/bin/bash');
        outScript.write("\n"+'cd '+CMSSWDir);
        outScript.write("\n"+'eval `scram runtime -sh`');
        outScript.write("\n"+'cd '+currentDir);
        outScript.write("\n"+command);
        outScript.close();
        os.system("chmod 777 "+currentDir+"/"+fn+".csh");
        #command2 = "bsub -q 1nd -cwd "+currentDir+" "+currentDir+"/"+fn+".sh";
        #os.system(command2);
        #print command2

for a in range(len(category)):
	for i in range(len(nameData)):
		fn = "Job/Condor_"+nameData[i]+"_"+category[a];
		outScript = open(fn,"w");
		outScript.write('universe = vanilla');
		outScript.write('\nExecutable = Job/Job_'+nameData[i]+'_'+category[a]+'.csh');
		outScript.write('\nRequirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000');
		outScript.write('\nShould_Transfer_Files = YES');
		outScript.write('\nWhenToTransferOutput = ON_EXIT');
		outScript.write('\nOutput		= condor_logs/'+nameData[i]+'_$(Cluster)_$(Process).stdout');
		outScript.write('\nError		= condor_logs/'+nameData[i]+'_$(Cluster)_$(Process).stderr');
		outScript.write('\nLog		= condor_logs/'+nameData[i]+'_$(Cluster)_$(Process).log');
		outScript.write('\nQueue 1');
		outScript.close();
		command2 = "condor_submit "+currentDir+"/"+fn;
		os.system(command2);
		print command2
