R_HOME := 		$(shell R RHOME)

## include headers and libraries for R 
RCPPFLAGS := 		$(shell $(R_HOME)/bin/R CMD config --cppflags)
RLDFLAGS := 		$(shell $(R_HOME)/bin/R CMD config --ldflags)

## include headers and libraries for Rcpp interface classes
RCPPINCL := 		$(shell echo 'Rcpp:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RCPPLIBS := 		$(shell echo 'Rcpp:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)
LIBS =-Wl,--no-as-needed -L./../build -lfastpca
EXE_INCLUDEPATH = -I./ -I./../include 
CFLAGS=-O3 -std=c++0x 
fastRPCA.so : 
	PKG_CPPFLAGS="$(RCPPFLAGS) $(RCPPINCL) $(CFLAGS) $(EXE_INCLUDEPATH)" PKG_LIBS="$(RLDFLAGS) $(RCPPLIBS) $(CFLAGS) $(LIBS)" R CMD SHLIB fastRPCA.cpp 
