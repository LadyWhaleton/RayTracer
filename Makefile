CXX       = g++
CXXFLAGS  = -g -O2 -Wall -c
LDFLAGS   = -g -O2 -Wall -L/usr/X11R6/lib -lGL -lglut
RM        = rm -f
EXEC      = ray_tracer

all: $(EXEC)

$(EXEC): ray_tracer.o main.o
	$(CXX) $(LDFLAGS) $^ -o $@

ray_tracer.o: ray_tracer.cpp ray_tracer.h
	$(CXX) $(CXXFLAGS) $< -o $@

main.o: main.cpp ray_tracer.h
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	$(RM) *~ *.bak *.o $(EXEC)
