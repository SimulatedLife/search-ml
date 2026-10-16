# hpc/ – SEARCH-ML job runner for the web portal

| Path | Purpose |
|---|---|
| `install_hpc.sh` | Installer: environment, scheduler detection, configuration, SSH key, self-test. |
| `bin/searchml-hpc` | Job manager the web portal calls over SSH (`ping`, `prepare`, `put`, `submit`, `status`, `get`, `cancel`, `list`, `cleanup`, `selftest`). |
| `bin/searchml-gateway` | Forced command for the web server's SSH key (no shell access). |
| `lib/searchml_hpc.py` | Implementation (Python standard library only). |
| `templates/slurm.sh`, `templates/pbs.sh` | Batch scripts generated for each job. |
| `searchml-hpc.conf.example` | Configuration template (the installer writes `searchml-hpc.conf`). |

Start with:

```bash
bash hpc/install_hpc.sh
```

Full guide: [../docs/HPC_SETUP.md](../docs/HPC_SETUP.md) · protocol: [../docs/ARCHITECTURE.md](../docs/ARCHITECTURE.md)
