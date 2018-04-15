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
$(TDIR)/forceQuads:	$(ODIR)/forceQuads.o $(LDIR)/libwsserver.a
	$(CXX) -g -o $(TDIR)/forceQuads $(ODIR)/forceQuads.o \
		-L$(LDIR) -lwsserver -legads -lpthread -lnlopt -lz $(RPATH) -lm

$(ODIR)/forceQuads.o:	$(SDIR)/forceQuads.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h $(IDIR)/wsserver.h
	$(CC) -g -c $(COPTS) $(DEFINE) -I$(IDIR) $(SDIR)/forceQuads.c -o $(ODIR)/forceQuads.o

clean:
	-rm $(ODIR)/forceQuads.o

cleanall:	clean
	-rm $(TDIR)/forceQuads
