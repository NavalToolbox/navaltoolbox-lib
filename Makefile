.PHONY: help install install-rust install-python lint lint-rust lint-python test test-rust test-python docs clean \
        docker-wheel docker-lint docker-test docker-docs docker-all

# Default target
help:
	@echo "┌─────────────────────────────────────────────────────────┐"
	@echo "│             NavalToolbox-lib Build System               │"
	@echo "└─────────────────────────────────────────────────────────┘"
	@echo ""
	@echo "Local targets:"
	@echo "  help            - Show this help message"
	@echo "  install         - Install everything (Rust + Python)"
	@echo "  install-rust    - Build Rust library"
	@echo "  install-python  - Build and install Python wheel"
	@echo "  lint            - Lint all code (Rust + Python)"
	@echo "  lint-rust       - Lint Rust code (cargo clippy)"
	@echo "  lint-python     - Lint Python code (flake8)"
	@echo "  test            - Run all tests (Rust + Python)"
	@echo "  test-rust       - Run Rust tests only"
	@echo "  test-python     - Run Python tests only"
	@echo "  docs            - Build Sphinx documentation"
	@echo "  docs-open       - Build and open documentation"
	@echo "  clean           - Remove build artifacts"
	@echo ""
	@echo "Docker targets:"
	@echo "  docker-wheel    - Build Python wheel using Docker"
	@echo "  docker-test     - Run all tests in Docker"
	@echo "  docker-docs     - Build documentation in Docker"
	@echo "  docker-all      - Build everything in Docker"

# ============================================================================
# Local targets
# ============================================================================

# Install everything
install: install-rust install-python
	@echo "✅ NavalToolbox installed successfully!"

# Build Rust library
install-rust:
	@echo "🦀 Building Rust library..."
	cd rust && cargo build --release
	@echo "✅ Rust library built!"

# Build and install Python wheel
install-python:
	@echo "🐍 Building Python wheel..."
	cd python && maturin develop --release
	@echo "✅ Python package installed!"

# Lint everything
lint: lint-rust lint-python
	@echo "✅ All linting complete!"

# Lint Rust code
lint-rust:
	@echo "🔍 Linting Rust code..."
	cd rust && cargo clippy -- -W clippy::all
	@echo "✅ Rust linting complete!"

# Lint Python code
lint-python:
	@echo "🔍 Linting Python code..."
	cd python && python3 -m flake8 . --count --statistics
	@echo "✅ Python linting complete!"

# Run all tests
test: test-rust test-python
	@echo "✅ All tests passed!"

# Run Rust tests only
test-rust:
	@echo "🦀 Running Rust tests..."
	cd rust && cargo test
	@echo "✅ Rust tests passed!"

# Run Python tests only
test-python:
	@echo "🐍 Running Python tests..."
	cd python && python -m pytest tests/ -v
	@echo "✅ Python tests passed!"

# Build documentation
docs:
	@echo "📚 Building Sphinx documentation..."
	cd docs && make html
	@echo "✅ Documentation built: docs/build/html/index.html"

# Build and open documentation
docs-open: docs
	@echo "🌐 Opening documentation..."
	open docs/build/html/index.html

# Clean build artifacts
clean:
	@echo "🧹 Cleaning build artifacts..."
	rm -rf rust/target/
	rm -rf python/.venv/
	rm -rf docs/build/
	rm -rf dist/
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name .pytest_cache -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete 2>/dev/null || true
	@echo "✅ Clean complete!"

# ============================================================================
# Docker targets
# ============================================================================

# Build Python wheel using Docker
docker-wheel:
	@echo "🐳 Building Python wheel in Docker..."
	docker build -t navaltoolbox-builder --target builder -f Dockerfile .
	docker run --rm -v $(PWD)/dist:/output navaltoolbox-builder
	@echo "✅ Wheel built in dist/ directory"

# Run tests in Docker
docker-test:
	@echo "🐳 Running tests in Docker..."
	docker build -t navaltoolbox-test --target test -f Dockerfile .
	docker run --rm navaltoolbox-test

# Build documentation in Docker
docker-docs:
	@echo "🐳 Building documentation in Docker..."
	docker build -t navaltoolbox-docs --target docs -f Dockerfile .
	docker run --rm -v $(PWD)/docs/build:/app/docs/build navaltoolbox-docs
	@echo "✅ Documentation built: docs/build/html/index.html"

# Build everything in Docker
docker-all: docker-wheel docker-test docker-docs
	@echo "✅ Docker build complete!"
