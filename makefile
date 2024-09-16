# 定义编译器  
CXX = g++  
  
# 定义编译选项  
CXXFLAGS = -mavx2 -fopenmp -O3 -march=native
  
# 定义链接器选项和库  
LDFLAGS = -lfftw3 -lm -fopenmp
  
# 定义源文件和目标文件  
SRC = main.cpp
OBJ = $(SRC:.cpp=.o)  
TARGET = main  
  
# 默认目标  
all: $(TARGET)  
  
# 链接目标  
$(TARGET): $(OBJ)  
	$(CXX) $(OBJ) -o $@ $(LDFLAGS)  
  
# 编译规则  
%.o: %.cpp  
	$(CXX) $(CXXFLAGS) -c $< -o $@  
  
# 清理编译生成的文件  
clean:  
	rm -f $(OBJ) $(TARGET)  
  
# 打印帮助信息  
help:  
	@echo "Usage:"  
	@echo "  make        - Compile and link the program."  
	@echo "  make clean  - Remove compiled files."  
	@echo "  make help   - Show this help message."


