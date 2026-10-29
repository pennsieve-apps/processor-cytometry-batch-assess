.PHONY: help test fixture fixture-demo run clean

SERVICE_NAME ?= processor-cytometry-batch-assess
PYTHON ?= python3
PREP ?= ../workflow-immune-cytometry/data/work/prep

.DEFAULT: help

help:
	@echo "Make Help for $(SERVICE_NAME)"
	@echo ""
	@echo "make test          - pytest: the synthetic cohort through the processor, both input layouts (no Docker)"
	@echo "make fixture       - write the synthetic cohort to data/fixture/{linked,flat}/"
	@echo "make fixture-demo  - subsample fcs-prep output of the demo data (PREP=$(PREP)) to 20k events/file,"
	@echo "                     with a planted +0.4 shift on CD3/CD4/CD8 in batch03, into data/fixture-demo/"
	@echo "make run           - build the image and run it in Docker over data/fixture/linked/, writing to data/output/"
	@echo "make clean         - remove data/output/"
	@echo ""
	@echo "PYTHON=$(PYTHON) (set PYTHON=.venv/bin/python after pip install -r requirements-dev.txt)"

test:
	$(PYTHON) -m pytest -q tests

fixture:
	rm -rf data/fixture
	$(PYTHON) -m tests.make_fixture data/fixture --layout both

fixture-demo:
	rm -rf data/fixture-demo
	$(PYTHON) -m tests.make_fixture data/fixture-demo --from $(PREP) --events 20000 --shift batch03:CD3,CD4,CD8:0.4

run: clean fixture
	mkdir -p data/output
	docker compose -f docker-compose.yml down --remove-orphans
	docker compose -f docker-compose.yml build
	docker compose -f docker-compose.yml up --exit-code-from processor
	docker compose -f docker-compose.yml down --remove-orphans
	@echo "--- outputs:"
	@find data/output -type f | sort

clean:
	rm -rf data/output
