# DO NOT MERGE THIS PART WITH MASTER

TARGET=nucops
LOCALSRC=nucops2.c

# The following code should handle all the headers and source files
# with the exception of the one containing "main" (specified above)

IDIR=include
ODIR=obj
BDIR=bin
SDIR=src

CC=gcc
CFLAGS= -lm -I$(IDIR) -fopenmp -O3

CFILES := $(wildcard $(SDIR)/*.c)
OFILES := $(patsubst $(SDIR)/%.c,$(ODIR)/%.o,$(CFILES))
IFILES := $(wildcard $(IDIR)/*.h)

$(ODIR)/%.o: $(SDIR)/%.c $(IFILES)
	@mkdir -p $(@D)
	$(info "Compiling" $<)
	$(CC) -c -o $@ $< $(CFLAGS)

$(BDIR)/$(TARGET): $(OFILES) $(LOCALSRC)
	@mkdir -p $(@D)
	$(CC) -o $@ $^ $(CFLAGS)

