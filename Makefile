# Caliper: High-Throughput Computational Metrology Platform
# Production Makefile with full verification and test targets.

CC ?= clang
CFLAGS ?= -O3 -Wall -Wextra -Werror -std=c99 -Isrc/c
LDLIBS ?= -lm
FUTHARK ?= futhark
PYTHON ?= python3

BUILD_DIR ?= _build
BIN_DIR ?= bin
LIB_DIR ?= lib

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    SHLIB_EXT := dylib
    SHLIB_FLAGS := -dynamiclib
else
    SHLIB_EXT := so
    SHLIB_FLAGS := -shared -fPIC
endif

FUTHARK_SRCS := $(wildcard src/futhark/*.fut)
FUTHARK_TESTS := $(wildcard tests/futhark/test_*.fut)

.PHONY: all build clean test test-futhark test-c test-python bench verify info help

all: build futhark-check

help:
	@echo "Caliper Build System Targets:"
	@echo "  make build         - Compile native C library and CLI binary"
	@echo "  make futhark-check - Type-check all Futhark kernels with --Werror"
	@echo "  make test          - Execute all Futhark, C, and Python test suites"
	@echo "  make test-futhark  - Compile and execute all Futhark kernel tests"
	@echo "  make test-c        - Compile and run C unit tests"
	@echo "  make test-python   - Run Python unit and integration test suite"
	@echo "  make bench         - Run metrology throughput benchmarks"
	@echo "  make verify        - Run forbidden dash checks and full test matrix"
	@echo "  make clean         - Remove compiled objects and build directories"

$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(BIN_DIR)
	@mkdir -p $(LIB_DIR)

build: $(BUILD_DIR)
	@echo "==> Compiling native C library and CLI executable..."
	$(CC) $(CFLAGS) $(SHLIB_FLAGS) -fPIC src/c/caliper_io.c src/c/caliper_bridge.c -o $(LIB_DIR)/libcaliper.$(SHLIB_EXT) $(LDLIBS)
	$(CC) $(CFLAGS) src/c/caliper_io.c src/c/caliper_bridge.c src/c/caliper_main.c -o $(BIN_DIR)/caliper-cli $(LDLIBS)
	@echo "==> Build successful: $(BIN_DIR)/caliper-cli, $(LIB_DIR)/libcaliper.$(SHLIB_EXT)"

futhark-check:
	@echo "==> Type-checking Futhark kernels..."
	@for f in $(FUTHARK_SRCS); do \
		echo "Checking $$f..."; \
		$(FUTHARK) check --Werror $$f || exit 1; \
	done
	@echo "==> All Futhark kernels passed strict type checks."

test-futhark: $(BUILD_DIR)
	@echo "==> Running Futhark test suites..."
	@for t in $(FUTHARK_TESTS); do \
		echo "Compiling and running $$t..."; \
		$(FUTHARK) c $$t -o $(BUILD_DIR)/$$(basename $$t .fut) || exit 1; \
		$(BUILD_DIR)/$$(basename $$t .fut) </dev/null > $(BUILD_DIR)/$$(basename $$t .fut).out || exit 1; \
		grep -q "true" $(BUILD_DIR)/$$(basename $$t .fut).out || (echo "FAIL: $$t did not output true" && exit 1); \
		echo "PASS: $$t"; \
	done
	@echo "==> All Futhark unit test suites PASSED."

test-c: $(BUILD_DIR)
	@echo "==> Running C unit tests..."
	$(CC) $(CFLAGS) src/c/caliper_io.c tests/c/test_c_io.c -o $(BUILD_DIR)/test_c_io $(LDLIBS)
	$(BUILD_DIR)/test_c_io
	$(CC) $(CFLAGS) src/c/caliper_io.c src/c/caliper_bridge.c tests/c/test_c_bridge.c -o $(BUILD_DIR)/test_c_bridge $(LDLIBS)
	$(BUILD_DIR)/test_c_bridge
	@echo "==> All C unit tests PASSED."

test-python:
	@echo "==> Running Python unit and integration tests..."
	$(PYTHON) -m unittest discover -s tests/python -v
	@echo "==> All Python test suites PASSED."

test: test-futhark test-c test-python
	@echo "======================================================"
	@echo "  ALL TEST SUITES PASSED SUCCESSFULLY (Futhark, C, Python)"
	@echo "======================================================"

bench:
	@echo "==> Running metrology throughput benchmarks..."
	$(PYTHON) bench/bench_pipeline.py --quick

verify:
	@echo "==> Running verification linters..."
	$(PYTHON) tools/check_no_em_dash.py
	@$(MAKE) futhark-check
	@$(MAKE) test
	@echo "==> Hard verification pass PASSED: Zero errors, zero warnings."

clean:
	@echo "==> Cleaning build artifacts..."
	@rm -rf $(BUILD_DIR) $(BIN_DIR) $(LIB_DIR) artifacts/
	@find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
	@find . -name "*.pyc" -delete 2>/dev/null || true
	@echo "==> Clean complete."
