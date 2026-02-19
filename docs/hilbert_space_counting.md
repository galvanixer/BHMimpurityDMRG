# Hilbert Space Counting

This note records how to count Hilbert-space dimensions for the TwoBoson model used in this repository.

## 1) Full unconstrained Hilbert space

For a chain of length $L$ with local cutoffs $n_{\max,a}$ and $n_{\max,b}$:

$$
d_{\mathrm{local}} = (n_{\max,a}+1)(n_{\max,b}+1),
$$
$$
D_{\mathrm{full}} = d_{\mathrm{local}}^{\,L}.
$$

This is the total number of basis states/eigenstates if no `(N_a, N_b)` sector is fixed.

## 2) Fixed quantum-number sector `(N_a, N_b)`

With `conserve_qns=true`, the Hamiltonian block dimension in sector $(N_a,N_b)$ is

$$
D_{\mathrm{sector}} = D_a D_b,
$$

where
$$
D_a = [x^{N_a}]\left(1+x+x^2+\cdots+x^{n_{\max,a}}\right)^L,
$$
$$
D_b = [y^{N_b}]\left(1+y+y^2+\cdots+y^{n_{\max,b}}\right)^L.
$$

$[x^N]$ means "coefficient of $x^N$".

### Closed forms

- If the site cap does not bind ($n_{\max}\ge N$):
$$
D(N,L,n_{\max})=\binom{N+L-1}{L-1}.
$$

- General bounded case (inclusion-exclusion):
$$
D(N,L,n_{\max})
=
\sum_{j=0}^{\left\lfloor N/(n_{\max}+1)\right\rfloor}
(-1)^j
\binom{L}{j}
\binom{N-j(n_{\max}+1)+L-1}{L-1}.
$$

Then use:
$$
D_a = D(N_a,L,n_{\max,a}),\quad
D_b = D(N_b,L,n_{\max,b}),\quad
D_{\mathrm{sector}}=D_aD_b.
$$

## 3) Worked example

Given:
$$
L=32,\quad n_{\max,a}=5,\quad n_{\max,b}=3,\quad N_a=32,\quad N_b=3.
$$

Counts:
$$
D_a = 562875591270069891,\qquad
D_b = 5984,
$$
$$
D_{\mathrm{sector}} = 3368247538160098227744.
$$

So the fixed-sector Hilbert-space dimension is approximately:
$$
D_{\mathrm{sector}} \approx 3.368\times 10^{21}.
$$

For comparison, the full unconstrained space is:
$$
D_{\mathrm{full}}
= 24^{32}
= 146811384664566452713597726037899455366168576.
$$
