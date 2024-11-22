C-P Reach
---------------

Cyber Physical System Reachability Analyzer.

# Installation

* Install [Poetry](https://python-poetry.org/docs/)

* Install CDD Lib
```bash
sudo apt install libcdd-dev
```

* Build Poetry Package

```bash
git clone git@github.com:cognipilot/c_p_reach
cd c_p_reach
poetry install
```

# Running Commands

First start a poetry shell in the poetry virtual environment.

```bash
poetry shell
```
Next, call the reachability tool.
```bash
c_p_reach
```

# Roadmap

## Working
* Multirotor 3D support
* Support Rover

## TODO
* Auto-read gains from CasADi model
* Take vehicle trajectory as input
