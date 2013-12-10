APP = $(PWD)/..

BIN = $(APP)/bin

AR = $(APP)/lib/base.a

LIB_PATH = -L /usr/local/lib/

LIBS = -pthread -lrt -lm -ldl -lcoinasl -lcoinblas -lcoinhsl -lipopt -lcoinlapack -lcoinmetis -lcoinmumps -lipoptamplinterface  

INCLUDE = -I $(APP)/include/ \
		  -I /usr/local/include/coin/

LDFLAGS = -shared -g -Wall -fPIC

CXXFLAGS =-g -Wall -fPIC

CC = g++

OBJS = problem.o \
       solution.o \
	   main.o

TGTS = $(BIN)/mf

.PHONY : all clean
all: $(TGTS)
	
.cc.o:
	$(CC) $(CXXFLAGS) -c $< -o $@ $(INCLUDE)

$(TGTS): $(OBJS)
	$(CC) $(CXXFLAGS) $(OBJS) -o $(TGTS) $(INCLUDE) $(AR) $(LIB_PATH) $(LIBS)

clean:
	rm $(OBJS) $(TGTS)
