### CDK2

# pf fetch structure=1FIN,1AQ1,1HCL,2CCH,2FVD --type cif --outdir data/cdk2
# pf fetch structure=CDK2_HUMAN --af --type=cif --outdir data/cdk2
# pf protein-extract data/cdk2/1FIN.cif --chains=A --remove-water --output data/cdk2/1FIN_a.cif
# pf protein-extract data/cdk2/1AQ1.cif --chains=A,B --remove-water --output data/cdk2/1AQ1_ab.cif
# pf protein-extract data/cdk2/1HCL.cif --chains=A --remove-water --output data/cdk2/1HCL_a.cif
# pf protein-extract data/cdk2/2FVD.cif --chains=A,B --remove-water --output data/cdk2/2FVD_ab.cif
# pf protein-extract data/cdk2/2CCH.cif --chains=A,G --remove-water --output data/cdk2/2CCH_ag.cif
# pf align-view --indir data/cdk2 P24941_AF 1FIN_a 1AQ1_ab 1HCL_a 2FVD_ab 2CCH_ag

### TYK2

# pf fetch structure=TYK2_HUMAN --af --type=cif --outdir data/tyk2
# pf fetch structure=6NZP,4OLI,5C03,6X8F,6AAM,4GJ2,7K7O --type cif --outdir data/tyk2
# pf protein-extract data/tyk2/4GJ2.cif --chains=A,B --remove-water --output data/tyk2/4GJ2_ab.cif
# pf protein-extract data/tyk2/5C03.cif --chains=A,C --remove-water --output data/tyk2/5C03_ac.cif
# pf protein-extract data/tyk2/6NZP.cif --chains=A,C --remove-water --output data/tyk2/6NZP_ac.cif
# pf protein-extract data/tyk2/6X8F.cif --chains=A,C --remove-water --output data/tyk2/6X8F_ac.cif
# pf protein-extract data/tyk2/7K7O.cif --chains=A,C --remove-water --output data/tyk2/7K7O_ac.cif
# pf protein-extract data/tyk2/4OLI.cif --chains=A,B --remove-water --output data/tyk2/4OLI_ab.cif
# pf protein-extract data/tyk2/6AAM.cif --chains=A,B --remove-water --output data/tyk2/6AAM_ab.cif
# pf align-view --indir data/tyk2 P29597_AF 6NZP_ac 4OLI_ab 5C03_ac 6X8F_ac 6AAM_ab 4GJ2_ab 7K7O_ac --method number

### BRAF

# pf fetch structure=BRAF_HUMAN --af --type=cif --outdir data/braf
# pf fetch structure=3OG7,4MNF,5C9C,4KSP,6PP9,7MFE --type cif --outdir data/braf
# pf protein-extract data/braf/3OG7.cif --chains=A,C --remove-water --output data/braf/3OG7_ac.cif
# pf protein-extract data/braf/4MNF.cif --chains=A,C --remove-water --output data/braf/4MNF_ac.cif
# pf protein-extract data/braf/5C9C.cif --chains=A,C --remove-water --output data/braf/5C9C_ac.cif
# pf protein-extract data/braf/4KSP.cif --chains=A,C --remove-water --output data/braf/4KSP_ac.cif
# pf protein-extract data/braf/6PP9.cif --chains=A,C --remove-water --output data/braf/6PP9_ac.cif
# pf protein-extract data/braf/7MFE.cif --chains=A --remove-water --output data/braf/7MFE_a.cif
# pf align-view --indir data/braf P15056_AF 3OG7_ac 4MNF_ac 5C9C_ac 4KSP_ac 6PP9_ac 7MFE_a

### CYP

# pf fetch structure=CP3A4_HUMAN,CP2C8_HUMAN,CP2B6_HUMAN,CP46A_HUMAN --af --type=cif --outdir data/cyp
# pf align-view --indir data/cyp P08684_AF P10632_AF P20813_AF Q9Y6A2_AF

# pf fetch structure=5VC0,3NXU,8SPD,8SO2,1PQ2,3IBD,4ENH --type cif --outdir data/cyp
# pf protein-extract data/cyp/5VC0.cif --chains=A,B,C --remove-water --output data/cyp/5VC0_abc.cif
# pf protein-extract data/cyp/3NXU.cif --chains=A,G,H --remove-water --output data/cyp/3NXU_abh.cif
# pf protein-extract data/cyp/8SPD.cif --chains=A,B,C --remove-water --output data/cyp/8SPD_abc.cif
# pf protein-extract data/cyp/8SO2.cif --chains=A,B,C,D,E,F,G --remove-water --output data/cyp/8SO2_abcdefg.cif
# pf protein-extract data/cyp/1PQ2.cif --chains=A,D --remove-water --output data/cyp/1PQ2_ad.cif
# pf protein-extract data/cyp/3IBD.cif --chains=A,B,C,D,E --remove-water --output data/cyp/3IBD_abcde.cif
# pf protein-extract data/cyp/4ENH.cif --chains=A,B,C --remove-water --output data/cyp/4ENH_abc.cif
# pf align-view --indir data/cyp P08684_AF 1PQ2_ad 3IBD_abcde 3NXU_abh 4ENH_abc 5VC0_abc 8SO2_abcdefg 8SPD_abc

### 3CLプロテアーゼ / Mpro

# pf fetch structure=R1AB_SARS2 --type=fasta --outdir data/mpro
# pf fetch structure=6LU7,6M2N,6ZRT,7BUY,7JQ2,7RVM,7RN1,7WOF,8FIW,8UPW --type=cif --outdir data/mpro
# pf protein-extract data/mpro/6LU7.cif --chains=A,B,C --remove-water --output data/mpro/6LU7_abc.cif
# pf protein-extract data/mpro/6M2N.cif --chains=A,E --remove-water --output data/mpro/6M2N_ae.cif
# pf protein-extract data/mpro/6ZRT.cif --chains=A,C --remove-water --output data/mpro/6ZRT_ac.cif
# pf protein-extract data/mpro/7BUY.cif --chains=A --remove-water --output data/mpro/7BUY_a.cif
# pf protein-extract data/mpro/7JQ2.cif --chains=A,B --remove-water --output data/mpro/7JQ2_ab.cif
# pf protein-extract data/mpro/7RVM.cif --chains=A,B --remove-water --output data/mpro/7RVM_ab.cif
# pf protein-extract data/mpro/7RN1.cif --chains=A,D --remove-water --output data/mpro/7RN1_ad.cif
# pf protein-extract data/mpro/7WOF.cif --chains=A,B --remove-water --output data/mpro/7WOF_ab.cif
# pf protein-extract data/mpro/8FIW.cif --chains=A,C --remove-water --output data/mpro/8FIW_ac.cif
# pf protein-extract data/mpro/8UPW.cif --chains=A,B --remove-water --output data/mpro/8UPW_ab.cif
# pf align-view --indir data/mpro P0DTD1_AF 6LU7_abc 6M2N_ae 6ZRT_ac 7BUY_a 7JQ2_ab 7RVM_ab 7RN1_ad 7WOF_ab 8FIW_ac 8UPW_ab
