# Finite elements method solver
Small rust util that solves heat transfer differential equation using finite element method. Equation is as following:
\[-(k(x)u'(x))' = 100x\]
\[u(2) = 0\]
\[u(0) + u'(0) = 20\]

# Usage 
To build program use:
```
cargo build --release
```

To solve equation use:
```
cargo run --release <N>
```
where N is the number of base functions. If not given, the default value is 10.
