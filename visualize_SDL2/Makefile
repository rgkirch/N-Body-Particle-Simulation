CC = g++

OSTYPE = $(shell uname -s)
ifeq ($(OSTYPE), Darwin)
  ### Old configurations for the old SDL for Mac. Change for SDL2 here ###
  INCLUDES = -I/Library/Frameworks/SDL.framework/Versions/A/Headers/
  LIBS = -lSDLmain -lSDL
  FRAMEWORKS = -framework SDL -framework OpenGL -framework Cocoa
else
  ### Configurations for SDL2 for Ubuntu -- Tested with Ubuntu 14.04 ###
  INCLUDES = -I/usr/local/include/SDL2
  LIBS = -lSDL2main -lSDL2 -lGL -lGLU
  FRAMEWORKS =
endif

visualize: visualize.cpp
	$(CC) -o $@ $< $(INCLUDES) $(FRAMEWORKS) $(LIBS)

clean:
	rm visualize
