CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS ?=
LDLIBS ?= -lflint -lmpfr -lgmp

BUILD_DIR := build
CERTIFICATE := $(BUILD_DIR)/certify

.PHONY: all certify clean

all: $(CERTIFICATE)

$(BUILD_DIR):
	mkdir -p $@

$(CERTIFICATE): certify.cc functions.cc methods.cc \
		functions.h methods.h | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) certify.cc functions.cc methods.cc -o $@ $(LDLIBS)

certify: $(CERTIFICATE)
	./$(CERTIFICATE)

clean:
	rm -rf $(BUILD_DIR)
