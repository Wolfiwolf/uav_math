BLD_DIR=build
SRC_DIR=src

INCS=$(shell find src -name "*.h")
SRCS=$(shell find src -name "*.c")
OBJS=$(SRCS:%=$(BLD_DIR)/%.o)

TARGET=uav_math_example

$(TARGET): $(OBJS)
	gcc -o $@ $^ -lm

$(BLD_DIR)/%.c.o: %.c $(INCS)
	@mkdir -p $(@D)
	gcc -c $< -o $@  

clean:
	rm -rf $(BLD_DIR)
	rm $(TARGET)

print:
	@echo $(OBJS)


