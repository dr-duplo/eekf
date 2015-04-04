# static library
SRC_LIB		:= eekf.c eekf_mat.c
TARGET_LIB	:= libeekf.a
OBJS_LIB	:= ${SRC_LIB:.c=.o}

# example program
SRC_EXAMPLE		:= examples/eekf_example.c
TARGET_EXAMPLE	:= examples/eekf_example
OBJS_EXAMPLE	:= ${SRC_EXAMPLE:.c=.o} $(TARGET_LIB)

# build params
BUILD_DIR		:= ./build
SRC_DIR			:= ./src
INCLUDE_DIRS	:= ./includes
VPATH			:= src

include toolchain_gcc.mk

.PHONY: clean

all: $(TARGET_LIB) $(TARGET_EXAMPLE)

# eekf archive
$(TARGET_LIB): $(OBJS_LIB) 
	@echo "[AR] archiving $@"
	@$(AR) $(BUILD_DIR)/$(TARGET_LIB) $(addprefix $(BUILD_DIR)/, $(OBJS_LIB))

# example program
$(TARGET_EXAMPLE): $(OBJS_EXAMPLE) 
	@echo "[LD] linking $@"
	@$(CC) -o $(BUILD_DIR)/$(TARGET_EXAMPLE) $(addprefix $(BUILD_DIR)/, $(OBJS_EXAMPLE)) $(LDFLAGS)

# compile rule
%.o: %.c
	@echo "[CC] compiling $@"
	@mkdir -p $(BUILD_DIR)/$(dir $@)
	@$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/$@ $<

# clean up rule
clean:
	@echo "[CLEAN] cleaning build files"
	@$(RM) -r $(BUILD_DIR)
