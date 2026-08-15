mk_fs () {
    name=$1
    shift
    for id in $* ; do
        fs="$fs data/${name}/${id}_main.cif"
    done
    echo $fs
}

### CDK2

# pf fetch structures=1FIN,1AQ1,1HCL,2CCH,2FVD --type cif --output data/cdk2
# pf fetch structure=CDK2_HUMAN --af --type=cif --output data/cdk2/af.cif
# pf protein-extract data/cdk2/1FIN.cif --chains=A --remove-water --output data/cdk2/1FIN_main.cif
# pf protein-extract data/cdk2/1AQ1.cif --chains=A,B --remove-water --output data/cdk2/1AQ1_main.cif
# pf protein-extract data/cdk2/1HCL.cif --chains=A --remove-water --output data/cdk2/1HCL_main.cif
# pf protein-extract data/cdk2/2FVD.cif --chains=A,B --remove-water --output data/cdk2/2FVD_main.cif
# pf protein-extract data/cdk2/2CCH.cif --chains=A,G --remove-water --output data/cdk2/2CCH_main.cif
# pf align-view data/cdk2/af.cif $(mk_fs cdk2 1FIN 1AQ1 1HCL 2FVD 2CCH)

### TYK2

# pf fetch structure=TYK2_HUMAN --af --type=cif --output data/tyk2/af.cif
# pf fetch structures=6NZP,4OLI,5C03,6X8F,6AAM,4GJ2,7K7O --type cif --output data/tyk2
# pf protein-extract data/tyk2/4GJ2.cif --chains=A,B --remove-water --output data/tyk2/4GJ2_main.cif
# pf protein-extract data/tyk2/5C03.cif --chains=A,C --remove-water --output data/tyk2/5C03_main.cif
# pf protein-extract data/tyk2/6NZP.cif --chains=A,C --remove-water --output data/tyk2/6NZP_main.cif
# pf protein-extract data/tyk2/6X8F.cif --chains=A,C --remove-water --output data/tyk2/6X8F_main.cif
# pf protein-extract data/tyk2/7K7O.cif --chains=A,C --remove-water --output data/tyk2/7K7O_main.cif
# pf protein-extract data/tyk2/4OLI.cif --chains=A,B --remove-water --output data/tyk2/4OLI_main.cif
# pf protein-extract data/tyk2/6AAM.cif --chains=A,B --remove-water --output data/tyk2/6AAM_main.cif
# pf align-view data/tyk2/af.cif $(mk_fs tyk2 6NZP 4OLI 5C03 6X8F 6AAM 4GJ2 7K7O) --method number

### BRAF

# pf fetch structure=BRAF_HUMAN --af --type=cif --output data/braf/af.cif
# pf fetch structures=3OG7,4MNF,5C9C,4KSP,6PP9,7MFE --type cif --output data/braf
# pf protein-extract data/braf/3OG7.cif --chains=A,C --remove-water --output data/braf/3OG7_main.cif
# pf protein-extract data/braf/4MNF.cif --chains=A,C --remove-water --output data/braf/4MNF_main.cif
# pf protein-extract data/braf/5C9C.cif --chains=A,C --remove-water --output data/braf/5C9C_main.cif
# pf protein-extract data/braf/4KSP.cif --chains=A,C --remove-water --output data/braf/4KSP_main.cif
# pf protein-extract data/braf/6PP9.cif --chains=A,C --remove-water --output data/braf/6PP9_main.cif
# pf protein-extract data/braf/7MFE.cif --chains=A --remove-water --output data/braf/7MFE_main.cif
# pf align-view data/braf/af.cif $(mk_fs braf 3OG7 4MNF 5C9C 4KSP 6PP9 7MFE)

### CTP

# pf fetch structure=CP3A4_HUMAN --af --type=cif --output data/cyp/cp3a4_af.cif
# pf fetch structure=CP2C8_HUMAN --af --type=cif --output data/cyp/cp2c8_af.cif
# pf fetch structure=CP2B6_HUMAN --af --type=cif --output data/cyp/cp2b6_af.cif
# pf fetch structure=CP46A_HUMAN --af --type=cif --output data/cyp/cp46a_af.cif
pf align-view data/cyp/cp3a4_af.cif data/cyp/cp2c8_af.cif data/cyp/cp2b6_af.cif data/cyp/cp46a_af.cif

# pf fetch structures=5VC0,3NXU,8SPD,8SO2,1PQ2,3IBD,4ENH --type cif --output data/cyp
# pf protein-extract data/cyp/5VC0.cif --chains=A,B,C --remove-water --output data/cyp/5VC0_main.cif
# pf protein-extract data/cyp/3NXU.cif --chains=A,G,H --remove-water --output data/cyp/3NXU_main.cif
# pf protein-extract data/cyp/8SPD.cif --chains=A,B,C --remove-water --output data/cyp/8SPD_main.cif
# pf protein-extract data/cyp/8SO2.cif --chains=A,B,C,D,E,F,G --remove-water --output data/cyp/8SO2_main.cif
# pf protein-extract data/cyp/1PQ2.cif --chains=A,D --remove-water --output data/cyp/1PQ2_main.cif
# pf protein-extract data/cyp/3IBD.cif --chains=A,B,C,D,E --remove-water --output data/cyp/3IBD_main.cif
# pf protein-extract data/cyp/4ENH.cif --chains=A,B,C --remove-water --output data/cyp/4ENH_main.cif

pf align-view data/cyp/cp3a4_af.cif $(mk_fs cyp 5VC0 3NXU 8SPD 8SO2 1PQ2 3IBD 4ENH)
