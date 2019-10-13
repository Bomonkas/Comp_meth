TARGET = main.out

HDRS = \
	   project/include

SRCS = \
	   project/src/main.cpp \
	   project/src/iterative.cpp \
	   project/src/slau.cpp

.PHONY: all clean

all: $(SRCS)
	$(CXX) -Wall -Wextra -Werror -I $(HDRS) -o $(TARGET) $(CXXFLAGS) $(SRCS)

clean:
	rm -rf $(TARGET)

