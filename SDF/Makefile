BASE = ../../rai
NAME   = $(shell basename `pwd`)
OUTPUT = lib$(NAME).so

OPTIM=fast
DEPEND = Core Algo Kin Gui Geo

SRCS = $(shell find . -maxdepth 1 -name '*.cpp' )
OBJS = $(SRCS:%.cpp=%.o)

include $(BASE)/build/generic.mk
