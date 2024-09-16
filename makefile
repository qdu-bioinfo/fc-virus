# 定义变量
CXX = g++
CXXFLAGS = -Wall -g
SRC_DIR = code
BIN_DIR = bin
TARGET = $(BIN_DIR)/FC-Virus

# 列出所有源文件
SRCS = $(SRC_DIR)/GeneralSet.cpp \
       $(SRC_DIR)/ReadUtility.cpp \
       $(SRC_DIR)/KmerUtility.cpp \
       $(SRC_DIR)/KmerHash.cpp \
       $(SRC_DIR)/HomoKmer.cpp \
       $(SRC_DIR)/Consensus.cpp \
       $(SRC_DIR)/FC-Virus.cpp

# 列出所有对象文件
OBJS = $(SRCS:.cpp=.o)

# 默认目标
all: $(TARGET)

# 规则来构建目标
$(TARGET): $(OBJS)
	@mkdir -p $(BIN_DIR) # 确保 bin 目录存在
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

# 规则来构建每个对象文件
%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 清理目标
clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
