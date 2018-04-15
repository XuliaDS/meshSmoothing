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
$(TDIR)/Vvalence:	$(ODIR)/valence.o $(LDIR)/libwsserver.a
	$(CXX) -g -o $(TDIR)/valence $(ODIR)/valence.o \
		-L$(LDIR) -lwsserver -legads -lpthread -lnlopt -lz $(RPATH) -lm

$(ODIR)/valence.o:	$(SDIR)/valence.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h $(IDIR)/wsserver.h
	$(CC) -g -c $(COPTS) $(DEFINE) -I$(IDIR) $(SDIR)/valence.c -o $(ODIR)/valence.o

clean:
	-rm $(ODIR)/valence.o

cleanall:	clean
	-rm $(TDIR)/valence
