# CAB-align latest version
CAB-align is a residue-residue contact area based protein structure alignment program.  
This is the latest version of CAB-align. Some bugs were fixed from the original version.
## Reference
Terashi, G. and Takeda-Shitaka, M., 2015. CAB-Align: A flexible protein structure alignment method based on the residue-residue contact area. PloS one, 10(10), p.e0141440.

https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0141440

## Set Up programs and Compile codes
Please compile codes in ```CAalign_src``` and ```CAMTX_src``` by  
```make``` command.

In ```cabalign.sh```, please specify the directory of the codes.  
example:  
```L9   BIN_DIR="/home/user1/CAB_ALIGN2"```

## Usage
```
cabalign.sh [PDB1 file] [PDB2 file] 
```

Output will be

```
...
No1 sco= 19650.751 rate= 0.580 CAD= 0.576 Sep>4 sco= 6421.931 rate= 0.501 CAD= 0.506 Weighted sco= 13036.341 rate= 0.558 CAD= 0.557 SupRec= 26.214 RMSD= 21.918 (269 Res)
-------------------------------------------------------------------------------------------DTILLTGLFAAFFTTFAFAPQSIKTIRTRNTEGISVVMYIMFLTGVISWIAYGIMRSDFAVLIANIVTLFLAAPVLVITLINRRKK------MDTILLTGLFAAFFTTFAFAPQSIKTIRTRNTEGISVVMYIMFLTGVISWIAYGIMRSDFAVLIANIVTLFLAAPVLVITLINRRKKHVLESSDTILLTGLFAAFFTTFAFAPQSIKTIRTRNTEGISVVMYIMFLTGVISWIAYGIMRSDFAVLIANIVTLFLAAPVLVITLINRRKKHVLESSG
DTILLTGLFAAFFTTFAFAPQSIKTIRTRNTEGISVVMYIMFLTGVISWIAYGIMRSDFAVLIANIVTLFLAAPVLVITLINRRKKHVLESDTILLTGLFAAFFTTFAFAPQSIKTIRTRNTEGISVVMYIMFLTGVISWIAYGIMRSDFAVLIANIVTLFLAAPVLVITLINRRKKHVLESSMDTILLTGLFAAFFTTFAFAPQSIKTIRTRNTEGISVVMYIMFLTGVISWIAYGIMRSDFAVLIANIVTLFLAAPVLVITLINRRKKHVLESMDTILLTGLFAAFFTTFAFAPQSIKTIRTRNTEGISVVMYIMFLTGVISWIAYGIMRSDFAVLIANIVTLFLAAPVLVITLINRRKKHVLE---
#FINISHED TOTAL TIME= 1.349745
```

- ```SupRec= 26.214``` is the normalized similarity score. It is normalized by the size of the internal residue-residue contact area and the contact area based similarity between the two input models.
