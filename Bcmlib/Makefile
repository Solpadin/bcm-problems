f       = f
c       = cpp
o       = o
a       = a
r       =

.SUFFIXES: .$f .$c .$o

CC      = mpicxx -O3 -c -I../qd
CF      = mpif77 -O3 -c 
LINK    =

LIB     = bcmlib.$a
LIBOBJS = \
	./Utils.$o ./Unit_mes.$o ./Kernel.$o ./Nurbs.$o \
	./Nodes.$o ./Cgrid.$o ./Cgrid_el.$o ./Cgrid_qg.$o ./Ccells.$o \
	./Shapes.$o ./Cacou.$o ./Cmapi.$o ./Cwave.$o \
	./Drafts.$o ./Ccohes2d.$o ./Ccohes3d.$o ./Cmindl2d.$o ./Cmindl3d.$o ./Clame2d.$o ./Clame3d.$o ./Cvisco2d.$o ./Cvisco2d_grad.$o\
	./Cacou3d.$o ./Cporosity2d.$o \
#	./CProfile.$o \
#	./CPro_bars.$o \
#	./CPro_bpr.$o \
#	./CPro_cable.$o \
#	./CPro_fitting.$o \
#	./CPro_lpr.$o \
#	./CPro_skin.$o \
#	./CPro_tower.$o \
#	./CPro_wire.$o \
#	./Montage.$o \
#	./Bсm_ident.$o \

$(LIB): $(LIBOBJS)
	ar r $(LIB) $(LIBOBJS)
	ranlib $(LIB)

.$f.$o:
	$(CF) $<

.$c.$o:
	$(CC) $<

clean:
	-rm *.o
	-rm $(LIB)


