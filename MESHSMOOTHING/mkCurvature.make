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
$(TDIR)/vCurvature:	$(ODIR)/vCurvature.o $(LDIR)/libwsserver.a
	$(CXX) -o $(TDIR)/vCurvature $(ODIR)/vCurvature.o \
		-L$(LDIR) -lwsserver -legads -lpthread -lz $(RPATH) -lm

$(ODIR)/vCurvature.o:	$(SDIR)/vCurvature.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h $(IDIR)/wsserver.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) $(SDIR)/vCurvature.c -o $(ODIR)/vCurvature.o

clean:
	-rm $(ODIR)/vCurvature.o

cleanall:	clean
	-rm $(TDIR)/vCurvature
