# Renamed cardiac *Pde classes to be *CellCollection
s/39 pack<void (MonodomainPde</42 pack<void (MonodomainTissue</
s/35 pack<void (BidomainPde</38 pack<void (BidomainTissue</
# Changed various cell models to be generated on-the-fly from CellML,
# which changed the class name
s/26 LuoRudyIModel1991OdeSystem/25 CellLuoRudy1991FromCellML/
s/30 BackwardEulerLuoRudyIModel1991/38 CellLuoRudy1991FromCellMLBackwardEuler/
s/33 BackwardEulerFoxModel2002Modified/39 CellFoxModel2002FromCellMLBackwardEuler/
s/29 DiFrancescoNoble1985OdeSystem/34 CellDiFrancescoNoble1985FromCellML/
s/20 Mahajan2008OdeSystem/25 CellMahajan2008FromCellML/
s/24 TenTusscher2006OdeSystem/32 CellTenTusscher2006EpiFromCellML/
s/21 Maleckar2009OdeSystem/26 CellMaleckar2008FromCellML/
s/44 HodgkinHuxleySquidAxon1952OriginalOdeSystem/31 CellHodgkinHuxley1952FromCellML/
s/21 FaberRudy2000Version3/27 CellFaberRudy2000FromCellML/
s/30 FaberRudy2000Version3Optimised/30 CellFaberRudy2000FromCellMLOpt/
# Changed a parameter name in the LuoRudy1991 cell model
s/31 fast_sodium_current_conductance/40 membrane_fast_sodium_current_conductance/
