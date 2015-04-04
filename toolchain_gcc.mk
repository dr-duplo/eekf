CC = gcc
AR = ar rcs
RM = rm -f

CFLAGS += -Wall -O2 -std=gnu99
CFLAGS += $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir))
LDFLAGS += -lm