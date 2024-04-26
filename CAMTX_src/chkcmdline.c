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
	cmd->DemoMode=false;
	cmd->ChkMissMode=false;
	cmd->r=1.2;
	cmd->MetaC=1.5;
	cmd->Color=0;
        p=argv;
        while (--argc){
         p++;
         if(**p == '-'){
          	switch(*(*p+1)){
		 case 'i':
			strcpy(cmd->p1,*(++p));
               		--argc; break;
		 case 'g':
			cmd->g_step=atof(*(++p));
               		--argc; break;
		 case 'm':
			cmd->MtxMode=true; 
			cmd->PlyMode=false; 
			//puts("#MidPointMode >> ON");
			break;
		 case 'p':
			cmd->PlyMode=true; 
			cmd->MtxMode=false; 
			//puts("#MidPointMode >> ON");
			break;
		 case 'P':
			cmd->ShowResNum=atoi(*(++p));
			cmd->DemoMode=true;
			cmd->PlyMode=false; 
			cmd->MtxMode=false; 
			//puts("#MidPointMode >> ON");
			--argc;
			break;
		 case 'k':
			cmd->ChkMissMode=true; 
			break;
		 case 'c':
			cmd->Color=atoi(*(++p));
			--argc;
			break;
		 case 'C':
			cmd->MetaC=atof(*(++p));
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
	puts("Usage: CAMTX -i [InputPDB file file] [(Options)]");
	puts("Marching cubes method");
	puts("Calculate Contact Area Matrix");
	puts("Options:");
	puts("-g [float] : Grid step size (def=1.0)");
	//puts("-r [float] : Radius of Solvent (def=1.2)");
	puts("-p 	 : Out put ply format file");
	puts("-P [int]   : Out put ply format for res_num=[int]");
	puts("-c [int]	 : coloring mode; 0: all white(=def), 1:residue, 2:contact");
	puts("-m 	 : Out put Contact Area Matrix(default=true)");
	puts("-C [float] : Metaball smoothing param (def=1.5)");
	puts("-k         : check missing residues (= no CA atom)");
	puts("Ver 0.20");
}
