hadd -f WWTree_data_golden_2p1.root WWTree_data_mu_*runD*
hadd -f WWTree_TTbar.root WWTree_TTbar_powheg_*
hadd -f WWTree_STop.root WWTree_*ch*root
hadd -f WWTree_WZ_excl.root WWTree_WZ_excl_*
hadd -f WWTree_ZZ_excl.root WWTree_ZZ_excl_*
mv WWTree_WJets.root WWTree_WJets_incl.root
hadd -f WWTree_WJets.root WWTree_WJets100.root WWTree_WJets200.root WWTree_WJets400.root WWTree_WJets600bis.root WWTree_WJets800.root WWTree_WJets1200.root WWTree_WJets2500.root 
hadd -f WWTree_VV.root WWTree_WW_excl.root WWTree_WZ_excl.root WWTree_ZZ_excl.root 
hadd -f WWTree_pseudodata.root WWTree_WJets.root WWTree_TTbar.root WWTree_VV.root WWTree_STop.root