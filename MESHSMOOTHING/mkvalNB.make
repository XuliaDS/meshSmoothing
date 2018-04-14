#
IDIR = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR = $(ESP_ROOT)/lib
ifdef ESP_BLOC
ODIR = obj
TDIR = testBUILD
else
ODIR = .
TDIR = $(ESP_ROOT)/bin
endif

SDIR = src

$(TDIR)/valenceNB:	$(ODIR)/valenceNB.o $(ODIR)/IO.o $(LDIR)/libwsserver.a
	$(CXX) -g -o $(TDIR)/valenceNB $(ODIR)/valenceNB.o $(ODIR)/IO.o\
		-L$(LDIR) -lwsserver -legads -lpthread -lnlopt -lz $(RPATH) -lm

$(ODIR)/IO.o:	$(SDIR)/IO.c $(IDIR)/egads.h $(SDIR)/IO.h $(IDIR)/egadsTypes.h \
	$(IDIR)/egadsErrors.h $(IDIR)/wsserver.h
		$(CC) -g -c $(COPTS) $(DEFINE) -I$(IDIR) $(SDIR)/IO.c   -o $(ODIR)/IO.o

$(ODIR)/valenceNB.o:	$(SDIR)/valenceNoBounds.c $(IDIR)/egads.h $(SDIR)/IO.h $(IDIR)/egadsTypes.h \
	$(IDIR)/egadsErrors.h $(IDIR)/wsserver.h
	$(CC) -g -c $(COPTS) $(DEFINE) -I$(IDIR)  $(SDIR)/valenceNoBounds.c -o $(ODIR)/valenceNB.o 

clean:
	-rm $(ODIR)/valenceNB.o $(ODIR)/IO.o

cleanall:	clean
	-rm $(TDIR)/valenceNB
