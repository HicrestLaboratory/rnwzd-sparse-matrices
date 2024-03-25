CC=c++
LIB_HOME=/path/to/lib ###########
LIBS= -pthread
INCLUDE=-I.
OPT= -std=c++20 -O3


MAIN=main.cpp *.cpp

########################

BUILDDIR:=build
TARGETDIR:=bin


all:$(TARGETDIR)/main


$(TARGETDIR)/main:$(MAIN) $(OBJECTS)
	@mkdir -p $(@D)
	$(CC) $^ -o $@ $(INCLUDE) $(LIBS) $(OPT)

clean:
	rm $(BUILDDIR)/* $(TARGETDIR)/* 

