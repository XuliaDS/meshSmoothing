#
IDIR = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR = $(ESP_ROOT)/lib
ifdef ESP_BLOC
ODIR = obj
TDIR = builds
else
ODIR = .
TDIR = $(ESP_ROOT)/bin
endif

SDIR = src
$(TDIR)/vFace:	$(ODIR)/vFace.o $(LDIR)/libwsserver.a
	$(CXX) -o $(TDIR)/vFace $(ODIR)/vFace.o \
		-L$(LDIR) -lwsserver -legads -lpthread -lz $(RPATH) -lm

$(ODIR)/vFace.o:	$(SDIR)/vFace.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h $(IDIR)/wsserver.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) $(SDIR)/vFace.c -o $(ODIR)/vFace.o

clean:
	-rm $(ODIR)/vFace.o

cleanall:	clean
	-rm $(TDIR)/vFace
