NAME := distributions-allhistos
SOURCE_JIT := $(NAME)-jit.cxx
SOURCE_NOJIT := $(NAME)-nojit.cxx

ROOT_CONFIG := `root-config --libs --cflags`

all: $(NAME)-jit $(NAME)-jit-opt $(NAME)-nojit $(NAME)-nojit-opt

$(NAME)-jit: $(SOURCE_JIT)
	g++ -O0 -o $(NAME)-jit $(SOURCE_JIT) $(ROOT_CONFIG) 

$(NAME)-jit-opt: $(SOURCE_JIT)
	g++ -O3 -o $(NAME)-jit-opt $(SOURCE_JIT) $(ROOT_CONFIG)

$(NAME)-nojit: $(SOURCE_NOJIT)
	g++ -O0 -o $(NAME)-nojit $(SOURCE_NOJIT) $(ROOT_CONFIG)

$(NAME)-nojit-opt: $(SOURCE_NOJIT)
	g++ -O3 -o $(NAME)-nojit-opt $(SOURCE_NOJIT) $(ROOT_CONFIG)

clean:
	rm -f *.so *.o $(NAME)-jit $(NAME)-jit-opt $(NAME)-nojit $(NAME)-nojit-opt
