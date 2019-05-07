#python python/produceWWNtuples.py -i eos/cms/store/user/arapyan/Run2 -n WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8  -o outputfile_EWK     -w 0.02696 -no 149993 -noNeg 0 -lumi 35900.0 --ismc 1 -trig 1 
#python python/produceWWNtuples.py -i eos/cms/store/user/arapyan/Run2 -n WpWpJJ_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8  -o outputfile_QCD -w 0.02615 -no 146429 -noNeg 0 -lumi 35900.0    --ismc 1 -trig 1

#python  python/produceWWNtuples.py -i /eos/cms/store/user/arapyan/mc  -n aqgc_ft1_2_allweights -o outputfile_aqgc_weighted_2p0     -w 0.02696 -no 1499993 -noNeg 0 -lumi 35900.0 --ismc 1 -trig 1
python  python/produceWWNtuples.py -i /eos/cms/store/user/arapyan/mc  -n aqgc_ft1_0p2_allweights -o outputfile_aqgc_unweighted_0p2_NoMwwCut   -w 0.02696 -no 1499993 -noNeg 0 -lumi 35900.0 --ismc 1 -trig 1

#python  python/produceWWNtuples.py -i /eos/cms/store/user/arapyan/mc  -n aqgc_fm1_5_allweights -o outputfile_aqgc_weighted_fm1_5p0     -w 0.02696 -no 1499993 -noNeg 0 -lumi 35900.0 --ismc 1 -trig 1

#python  python/produceWWNtuples.py -i /eos/cms/store/user/arapyan/mc  -n aqgc_fs0_5_allweights -o outputfile_aqgc_weighted_fs0_5p0     -w 0.02696 -no 1499993 -noNeg 0 -lumi 35900.0 --ismc 1 -trig 1
#//WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8
#//WWJJ_SS_WToLNu_EWK-QCD_aQGC-FT-FS-FM_TuneCUETP8M1_13TeV_madgraph-pythia8
#//WpWpJJ_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8

