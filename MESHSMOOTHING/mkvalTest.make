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
$(TDIR)/valence:	$(ODIR)/valenceTest.o $(LDIR)/libwsserver.a
	$(CXX) -g -o $(TDIR)/valenceTest $(ODIR)/valenceTest.o \
		-L$(LDIR) -lwsserver -legads -lpthread -lnlopt -lz $(RPATH) -lm

$(ODIR)/valenceTest.o:	$(SDIR)/valenceTest.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h $(IDIR)/wsserver.h
	$(CC) -g -c $(COPTS) $(DEFINE) -I$(IDIR) $(SDIR)/valenceTest.c -o $(ODIR)/valenceTest.o

clean:
	-rm $(ODIR)/valenceTest.o

cleanall:	clean
	-rm $(TDIR)/valenceTest
