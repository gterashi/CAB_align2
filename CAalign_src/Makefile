#OBJECTS = chkcmdline.o readMAT.o readpdb.o sch.o Gcounter.o TOZ.o KAITEN.o VEC_kei.o kaiten.o cd_func.o SF_counter.o sch_tri.o Mcounter.o MVandSCORE.o lsfit.o det_sf.o enephobic_ca.o high_reso_doc.o

CC = gcc
#COM =  readpdb.o  cd_func.o chkcmdline.o rmsd.o
COM =  chkcmdline.o readpdb.o cd_func.o dp.o rmsd.o mammoth.o
OBJECTS = dali.o tm.o $(COM)
OBJECTS2 = shib_func_table.o $(COM)
HEADER = struct.h func.h dali.h tm.h
PG = CAalign
PG2 = shib_ANA
PG2 = shib_A5_simple

#OPT = -O3 -msse -msse2 -msse3 -mmmx -ftree-vectorize -ftree-vectorizer-verbose=3 -ffast-math -march=native -funroll-loops -fopenmp
OPT = -O3 -msse -msse2 -msse3 -mmmx -ftree-vectorize -ftree-vectorizer-verbose=3 -ffast-math  -funroll-loops 

$(PG) : main.c $(OBJECTS) $(HEADER) $(COM)
	$(CC) -o $(PG) main.c $(OBJECTS) -lm $(OPT) 

$(PG2) : main.c $(OBJECTS2) $(HEADER) $(COM)
	$(CC) -o $(PG2) main.c $(OBJECTS2) -lm $(OPT)

dali.o : dali.c $(HEADER) 
	$(CC) -c dali.c $(OPT)

tm.o : tm.c $(HEADER) -lm
	$(CC) -c tm.c $(OPT)

mammoth.o : mammoth.c $(HEADER) -lm
	$(CC) -c mammoth.c $(OPT)

pdbfft.o : pdbfft.c $(HEADER)
	$(CC) -c pdbfft.c $(OPT)

rmsd.o : rmsd.c rmsd.h
	$(CC) -c rmsd.c $(OPT) -lm

chkcmdline.o : chkcmdline.c $(HEADER)
	$(CC) -O3 -c chkcmdline.c

shib_func.o : shib_func.c $(HEADER)
	$(CC) -O3 -c shib_func.c -lm $(OPT)

shib_func_table.o : shib_func_table.c $(HEADER)
	$(CC) -O3 -c shib_func_table.c -lm $(OPT)

dp.o : dp.c $(HEADER)
	$(CC) -O3 -c dp.c -lm $(OPT)

chkcmdline_clst.o : chkcmdline_clst.c $(HEADER)
	$(CC) -O3 -c chkcmdline_clst.c

side.o : side.c $(HEADER)
	$(CC) -O3 -c side.c

ENE.o : ENE.c $(HEADER)
	$(CC) -O3 -c ENE.c

cd_func.o : cd_func.c $(HEADER)
	$(CC) $(OPT) -c cd_func.c

mk_trunk_rescd.o : mk_trunk_rescd.c $(HEADER)
	$(CC) -O3 -c mk_trunk_rescd.c

chkcmdline_st.o : chkcmdline_st.c $(HEADER)
	$(CC) -O3 -c chkcmdline_st.c

CHK.o : CHK.c $(HEADER)
	$(CC) -O3 -c CHK.c

hbond.o : hbond.c $(HEADER)
	$(CC) -O3 -c hbond.c

SA.o : SA.c $(HEADER)
	$(CC) -O3 -c SA.c

readMESH.o : readMESH.c $(HEADER)
	$(CC) -O3 -c readMESH.c

readpdb.o : readpdb.c $(HEADER)
	$(CC) -O3 -c readpdb.c

lsfit.o : lsfit.c $(HEADER)
	$(CC) -O3 -c lsfit.c $(OPT)

maxsub.o : maxsub.c $(HEADER)
	$(CC) -O3 -c maxsub.c

cluster.o : cluster.c $(HEADER)
	$(CC) -O3 -c cluster.c

maxsub_log.o : maxsub_log.c $(HEADER)
	$(CC) -O3 -c maxsub_log.c

get_chi.o : get_chi.c $(HEADER)
	$(CC) -O3 -c get_chi.c
