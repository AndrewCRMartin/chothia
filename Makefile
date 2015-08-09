COPT	= -Wall -ansi -I$(HOME)/include
LINK1	= -L$(HOME)/lib -lbiop -lgen -lxml2
LINK2	=
CC	= cc

EXE	= chothia
OFILES	= chothia.o KabCho.o
LFILES  = 

$(EXE) : $(OFILES) $(LFILES)
	$(CC) -o $(EXE) $(OFILES) $(LFILES) $(LINK1) -lm $(LINK2)

.c.o :
	$(CC) $(COPT) -o $@ -c $<

clean :
	/bin/rm -f $(EXE) $(OFILES) $(LFILES)
