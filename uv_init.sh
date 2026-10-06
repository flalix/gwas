# 1. Pin the interpreter for this project (writes .python-version)
uv python install 3.10        # downloads a managed CPython if you don't have one
uv python pin 3.10  # chose the correct

# 2. Create the venv and resolve/install everything
cd /home/flavio/uv/gwas

uv sync --locked

.venv/bin/python -c "import sys, ipykernel; print(sys.executable); print(sys.version)"

.venv/bin/python -m ipykernel install \
  --user \
  --name gwas \
  --display-name "GWAS (Python 3.10)"

# 3. Sanity check
uv run python -c "import yaml, numpy, dotenv, kaleido; print('ok')"
