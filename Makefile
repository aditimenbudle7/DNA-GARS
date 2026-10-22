# DNA-GARS build
CC      := gcc
CFLAGS  := -std=gnu11 -Wall -Wextra -O2 -Iinclude
LDFLAGS :=

SRC_DIRS := src src/search src/align src/reconstruct
SRCS     := $(foreach d,$(SRC_DIRS),$(wildcard $(d)/*.c))
OBJS     := $(patsubst %.c,build/%.o,$(SRCS))
LIB_OBJS := $(filter-out build/src/main.o,$(OBJS))

TARGET := dna-gars
TESTS  := test_search test_align test_reconstruct

.PHONY: all test clean debug
.SECONDARY:

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

build/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

# Each test links its own main with every module except src/main.c
test_%: build/tests/test_%.o $(LIB_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

test: $(TESTS)
	@for t in $(TESTS); do echo "== $$t =="; ./$$t || exit 1; done

# Debug build with AddressSanitizer + UBSan (catches memory bugs)
debug: CFLAGS  := -std=gnu11 -Wall -Wextra -g -O0 -fsanitize=address,undefined -Iinclude
debug: LDFLAGS += -fsanitize=address,undefined
debug: clean $(TARGET)

clean:
	rm -rf build $(TARGET) $(TESTS)

-include $(OBJS:.o=.d)