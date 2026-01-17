# Multi-stage Dockerfile for navaltoolbox-lib
# Supports: wheel building, testing, and documentation
#
# Stages:
#   - rust-builder: Builds the Rust library
#   - builder: Builds the Python wheel with maturin
#   - test: Runs all tests (Rust + Python)
#   - docs: Builds Sphinx documentation
#   - final: Outputs the wheel package

# ============================================================================
# Base Rust stage
# ============================================================================
FROM rust:1.75-slim as rust-base

RUN apt-get update && apt-get install -y \
    build-essential \
    pkg-config \
    python3 \
    python3-pip \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy Rust project files
COPY rust/ ./rust/
COPY python/ ./python/
COPY README.md ./

# ============================================================================
# Rust builder stage - compiles Rust library
# ============================================================================
FROM rust-base as rust-builder

RUN echo "🦀 Building Rust library..." && \
    cd rust && cargo build --release && \
    echo "✅ Rust build complete!"

# ============================================================================
# Python wheel builder stage - builds wheel with maturin
# ============================================================================
FROM rust-base as builder

# Install maturin
RUN pip3 install --break-system-packages maturin

RUN echo "📦 Building Python wheel..." && \
    cd python && maturin build --release && \
    echo "✅ Wheel built:" && \
    ls -lh ../rust/target/wheels/

# Default: copy wheel to output
WORKDIR /output
CMD ["sh", "-c", "cp /app/rust/target/wheels/*.whl /output/ && ls -lh /output/"]

# ============================================================================
# Test stage - runs all tests
# ============================================================================
FROM rust-base as test

# Install test dependencies
RUN pip3 install --break-system-packages maturin pytest

# Build and install
RUN cd python && maturin develop --release

# Copy test files
COPY rust/tests/ ./rust/tests/

# Run Rust tests
RUN echo "🦀 Running Rust tests..." && \
    cd rust && cargo test && \
    echo "✅ Rust tests passed!"

# Run Python tests
RUN echo "🐍 Running Python tests..." && \
    cd python && python3 -m pytest tests/ -v && \
    echo "✅ Python tests passed!"

CMD ["sh", "-c", "cd rust && cargo test && cd ../python && python3 -m pytest tests/ -v"]

# ============================================================================
# Documentation stage - builds Sphinx documentation
# ============================================================================
FROM python:3.11-slim as docs-base

RUN apt-get update && apt-get install -y \
    make \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    sphinx \
    pydata-sphinx-theme \
    sphinx-copybutton \
    sphinx-design

WORKDIR /app

# Copy documentation files
COPY docs/ ./docs/

FROM docs-base as docs

RUN echo "📚 Building Sphinx documentation..." && \
    cd docs && make html && \
    echo "✅ Documentation built!"

CMD ["sh", "-c", "cd docs && make html"]

# ============================================================================
# Final stage - outputs wheel package
# ============================================================================
FROM python:3.11-slim as final

WORKDIR /output

COPY --from=builder /app/rust/target/wheels/*.whl /dist/

CMD ["sh", "-c", "cp /dist/*.whl /output/ && ls -lh /output/"]
