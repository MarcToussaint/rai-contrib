BASE = ../../rai
NAME   = $(shell basename `pwd`)
OUTPUT = lib$(NAME).so

DEPEND = Core Algo Kin Gui Geo KOMO Optim

DEPEND_UBUNTU += libompl-dev
LIBS += -lompl

SRCS = $(shell find . -maxdepth 1 -name '*.cpp' )
OBJS = $(SRCS:%.cpp=%.o)

include $(BASE)/makeutils/generic.mk
