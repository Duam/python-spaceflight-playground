# Spaceflight-playground

This repository is a collection of spaceflight-related optimal-control problems
and visualizations.

#### Folder structure
```
./
 ├ apps/                          # Executable files
 ├ legacy/                        # Old matlab code
 ├ src/
 │ └ spaceflight_playground/      # Source code, library
 │   ├ orbiter/                       # Point-mass object in orbit
 │   ├ booster/                       # Rigid-body object during liftoff
 │   └ aux_models/                    # Auxiliary models
 └ tests/                         # Unit-tests
```

#### Usage
```bash
$ git clone git@github.com:Duam/python-spaceflight-playground.git
$ cd python-spaceflight-playground
$ virtualenv venv
$ source venv/bin/activate
$ pip install -e .
$ echo "Enjoy!"
```
