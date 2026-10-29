FROM python:3.12.11-slim-bookworm@sha256:519591d6871b7bc437060736b9f7456b8731f1499a57e22e6c285135ae657bf7

ENV PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /processor

# Pinned Python dependencies (wheels, and FlowIO's pure-Python sdist), plus the
# Lambda runtime interface client for running as a Lambda function.
COPY requirements.txt /processor/requirements.txt
RUN pip install -r /processor/requirements.txt awslambdaric==2.2.1

COPY processor/ /processor/processor/
COPY VERSION /processor/VERSION

ENV PYTHONPATH="/processor" \
    PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/matplotlib

# Build matplotlib's font cache now (as the run-time user would), so a run needs no writes outside /tmp.
RUN python -c "import matplotlib.pyplot" && chmod -R a+rwX /tmp/matplotlib

# Runs as uid 1000 on the platform; nothing here needs root.
USER 1000

# Runs main() once as a container, or serves processor.handler as a Lambda
# function (see processor/entry.py).
ENTRYPOINT ["python3", "-m", "processor.entry"]
