#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "struct.h"

int chkcmdline( int argc, char **argv,CMD *cmd){
        char **p;
	int num=0;
	void errmsg();
        if(argc < 3){
		errmsg();
                return(FALSE);
        }
	//default status
	cmd->g_step=1.00;
	cmd->MtxMode=true;
	cmd->PlyMode=false;
	cmd->IgNeighbor=true;
	cmd->CrossOver=false;
	cmd->r=1.2;
	cmd->z=1.0;
	cmd->zmc=1.0;
	cmd->MetaC=1.5;
	cmd->len=6;
	cmd->area=5.00;
	cmd->LocalRate=0.0;//Local structure similarity from DALI elastic score
	cmd->ExtRate=0.1;
	cmd->GapOpen=90.0;
	cmd->GapExt=0.0;
	cmd->SepWeight=0.5;
	cmd->LQmode=false;
	cmd->FilterOut=-1.00;
        p=argv;
        while (--argc){
         p++;
         if(**p == '-'){
          	switch(*(*p+1)){
		 case 'i':
			strcpy(cmd->d1,*(++p));
               		--argc; break;
		 case 'I':
			strcpy(cmd->d2,*(++p));
               		--argc; break;
		 case 'p':
			strcpy(cmd->p1,*(++p));
               		--argc; break;
		 case 'P':
			strcpy(cmd->p2,*(++p));
               		--argc; break;
		 case 'g':
			cmd->IgNeighbor=true;
               		break;
		 case 'u':
			cmd->IgNeighbor=false;
               		break;
		 case 'q':
			cmd->LQmode=true;
               		break;
		 case 'm':
			cmd->MtxMode=true; 
			cmd->PlyMode=false; 
			//puts("#MidPointMode >> ON");
			break;
		 case 'v':
			cmd->CrossOver=true; 
			break;
		 case 'C':
			cmd->MetaC=atof(*(++p));
               		--argc; break;
		 case 'L':
			cmd->len=atoi(*(++p));
               		--argc; break;
		 case 'z':
			cmd->z=atof(*(++p));
               		--argc; break;
		 case 'Z':
			cmd->zmc=atof(*(++p));
               		--argc; break;
		 case 'a':
			cmd->area=atof(*(++p));
               		--argc; break;
		 case 'l':
			cmd->LocalRate=atof(*(++p));
               		--argc; break;
		 case 'e':
			cmd->ExtRate=atof(*(++p));
               		--argc; break;
		 case 'G':
			cmd->GapOpen=atof(*(++p));
               		--argc; break;
		 case 'E':
			cmd->GapExt=atof(*(++p));
               		--argc; break;
		 case 's':
			cmd->SepWeight=atof(*(++p));
               		--argc; break;
		 case 'f':
			cmd->FilterOut=atof(*(++p));
               		--argc; break;
		 default: 
		  	fprintf(stderr,"No such a option: %s\n",*p+1); 
			errmsg(); 
			return(FALSE); 
		 break;
	 }
	 }
        }
	//option check
	 
	//printf("#File Name [%s]\n",cmd->p1);
	//printf("#Grid Step [%f]\n",cmd->g_step);
	//printf("#MetaballC [%f]\n",cmd->MetaC);
	//printf("#Sol Radius[%f]\n",cmd->r);
        return(TRUE);
}

void errmsg(){
	puts("Usage: CAalign -i [contact matrix1] -I [contact matrix2] -p [PDB1] -P [PDB2] [(Options)]");
	puts("Residue-residue contact area based Alignment");
	puts("Local Initial Alignments(mammoth, TMalign) -> Contact area based iterative DP");
	puts("Options:");
	puts("-g       : ignore +-1 neighbor residue contact based on real number (default=true)");
	puts("-u       : use +-1 neighbor residue contact based on real number (default=false)");
	//puts("-L int   : Length of a fragment (default= 6)");
	//puts("-a float : Cut off of Local Contact Area per residue(default=5)");
	//puts("-l float : Cut off of local structure similarity (default=0.0)");
	//puts("-e float : Cut off of extended CAD-score (default=0.1)");
	//puts("-z float : Z-score cutoff of doublet segments (default= 1.0)");
	//puts("-Z float : Z-score cutoff of Iteration (default= 1.0)");
	//puts("-v       : Allow cross Over mode, (sequence order independent)");
	puts("-G float : Gap open penalty in iter DP (default= 90.0)");
	puts("-E float : Gap extension penalty in iter DP (default= 0.0)");
	puts("-s float : Weight of near pair (i-j<5) (default= 0.5)");
	puts("-q       : Show Local Quality (default=false)");
	//puts("-f float : Filter Out Mode, ignore Max{SupR} < float");
	puts("***Ver 1.00***");
/*
	puts("Including TMalign algorith : 20131224");
	puts("Add Suboptimal Alignment Algorithm: 20140123");
	puts("refine Suboptimal Alignment Algorithm: 20140218");
	puts("Add sequence separation data(sep>=4) : 20140227");
	puts("Add sequence separation weight to scoring function: 20140304 v0.65");
	puts("Add Superfamily recognition score: 20140318 v0.66");
	puts("Add Local Quality score mode (-q): 20140410 v0.661");
	puts("Set optimal param to default param: 20141219 v1.00");
*/
}
