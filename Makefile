DIR      = .
TARGET   = "${DIR}/test"

CC       = /usr/csshare/pkgs/gcc-4.8.1/bin/g++
WFLAGS   = -W -Wall -Wextra -Werror -Wfatal-errors
CFLAGS   = -c -O2 -ansi -std=c++0x -I/home/rlair/include
LDFLAGS  = -L/home/rlair/lib
LDLIBS	= -lboost_serialization

HDRS := $(filter-out aki%, $(wildcard *.h))
SRCS := $(filter-out aki%, $(wildcard *.cpp))
OBJS = $(SRCS:.cpp=.o)

.PHONY: all clean build info rebuild

all: clean build

$(TARGET): $(OBJS)
	$(CC) $^ $(LDFLAGS) $(LDLIBS) -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) $<

info:
	$(CC) $(CFLAGS) -dM -E - < /dev/null

clean:
	rm -f $(OBJS) $(TARGET)

build: $(TARGET)

rebuild: clean build
