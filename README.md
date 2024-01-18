# Finite elements method solver
Small rust util that solves heat transfer differential equation using finite element method. Equation is as following: <br>
$-(k(x)u'(x))' = 100x$ <br>
$u(2) = 0$ <br>
$u(0) + u'(0) = 20$ <br>
$$k(x) =
 \begin{cases}
    1  & \text{dla} x\in [0,1]\\
    2x & \text{dla} x\in (1,2]
  \end{cases}$$

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

