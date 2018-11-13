## definitely to-do

- add tests
- make sure package passes `devtools::check()` etc.
- sort out issues with solvers: `rk4` vs. `lsoda` vs `ode45` vs ...

## API questions

- make link functions part of the model object, rather than the fit?
- difference operator as part of formula (`.diff()`) rather than additional argument?
- change name of `wmvrnorm`/synonymize to `impsamp`?

## probably doable

- incorporate `dnorm_n` (sd-profiled `dnorm`, from development version of `bbmle`) in examples
- priors/regularization/MAP estimation
- allow `predict` to use new predictor variables (e.g. finer time steps) [potentially tricky if diff-variables are being used ...]
- discussion of similar approaches/packages (nlmeODE, Stan + ODE module, odeintr, PKPD models ...)	

## ambitious

- method for using an ensemble with hypercube/Sobol starting values
    - plus plotting methods for Raue et al multi-mode detection (i.e., cumulative distribution of negative log-likelihoods)
- self-starting models (see `?selfStart`)
- vector-valued	parameters, state variables (e.g. for setting up complex compartmental models with the linear chain trick)

## almost certainly too hard

- allow gradient as well as trajectory matching
- C++/compiled version (e.g. via NIMBLE?)

## who knows?

- try out `nloptr::auglag` for profile CIs of predictions etc.?
- allow `predict` to use new parameter values (maybe not necessary since we have the `wvmrnorm` method, which is the main reason to want to do this)
- see whether the current importance-sampling machinery in the devel version of `bbmle` can replace `wmvrnorm` (or vice-versa) ...
- prettier plotting methods (with `ggplot2`)
