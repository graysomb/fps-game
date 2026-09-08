# Greedy cube activation demos

All three animations use the same 1,692-cell floating irregular structure and
600 fixed simulation steps. The first 20 steps show the static source.

- `greedy-cpu-st.gif`: greedy groups on the single-threaded CPU solver.
- `greedy-cpu-mt.gif`: the same groups on the worker-pool CPU solver.
- `fine-cpu-st.gif`: ordinary fine degrees of freedom with fracture disabled.

The greedy cover produced 113 groups with edge sizes 1 through 5. Its 2,324
shared child particles became 177 independent controls and 2,147 dependent
samples. CPU ST and MT matched exactly over 1,350,244 recorded states.
