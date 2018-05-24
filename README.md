# Dependencies
python3 with: casadi, numpy, matplotlib

# ToDo
- [ ] Convert existing matlab code to python
- - [x] Integrators
- - [x] Polar orbit model & test
- - [x] Cartesian orbit model & test
- - [ ] OCPs
- - [ ] Utils
- - [ ] Get data from old .mat files for comparing results
- [x] Test for liftoff model
- [ ] Check liftoff model equations (particularly torque and thrust)
- [x] Test for log wind model
- [x] Test for kepler orbit model
- [x] ~~Edge cases for kepler orbit model: Numerical errors at e -> 0 or e -> 1~~ Used a different equation for discretization
- [ ] Add a data folder for precomputed initial guesses (.csv)