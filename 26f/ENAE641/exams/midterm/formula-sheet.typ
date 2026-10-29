/////////////
// Imports //
/////////////
#import "@preview/physica:0.9.8": *
#import "@preview/unify:0.8.1": *

/////////////////
// Maths Setup //
/////////////////

// upright vectors
#let vectorboldupright(a) = vb($upright(#a)$)
#let vbu = vectorboldupright
#let vectorunitupright(a) = vu($upright(#a)$)
#let vuu = vectorunitupright
#let vectorarrowupright(a) = va($upright(#a)$)
#let vau = vectorarrowupright

// automatically use square brackets for vectors and matrices
#set math.vec(delim: "[")
#set math.mat(delim: "[")
#let vecrowOld = vecrow
#let vecrow = vecrowOld.with(delim: "[")

// custom units
#add-unit("decibelIsotropic", "dBi", "upright(\"dBi\")", space: true)
#add-unit("decibelWatt", "dBW", "upright(\"dBW\")", space: true)

// extra operators
#let tr = math.op("tr")
#let adj = math.op("Adj")
#let blkdiag = math.op("blkdiag")
#let spn = math.op("span")

//////////////////
// Sheet Layout //
//////////////////

// one US-letter sheet (two sides), three columns, tight margins
#set page(paper: "us-letter", margin: 0.25in, columns: 3)
#set columns(gutter: 7pt)
#set text(font: "New Computer Modern", size: 7.3pt)
#set par(leading: 0.35em, spacing: 0.45em)
#set list(indent: 0pt, body-indent: 0.3em, spacing: 0.35em)
#set enum(indent: 0pt, body-indent: 0.3em, spacing: 0.35em)
#set math.mat(row-gap: 0.2em, column-gap: 0.35em)
#show math.equation.where(block: true): set block(spacing: 0.35em)
#set table(inset: (x: 1.5pt, y: 1.2pt), stroke: 0.3pt + luma(170))
#set image(width: 100%)

// section colours (one hue per topic group)
#let c-la = rgb("#1c5cab")
#let c-mod = rgb("#006300")
#let c-sol = rgb("#4a3aa7")
#let c-exp = rgb("#a8440f")
#let c-eig = rgb("#b42828")
#let c-stab = rgb("#0f6e6e")

// section header bar
#let sec(title, c) = block(
    width: 100%,
    fill: c,
    inset: (x: 2pt, y: 1.8pt),
    above: 4.5pt,
    below: 3pt,
    sticky: true,
    text(fill: white, weight: "bold", size: 7.5pt, title),
)
// run-in sub-heading in the section colour
#let sub(c, body) = text(fill: c, weight: "bold", body)
// boxed must-know result
#let key(body) = box(stroke: 0.5pt + luma(80), inset: (x: 1.5pt, y: 1.5pt), outset: (y: 0.5pt), body)
// "careful / not true in general"
#let warn(body) = block(
    width: 100%,
    fill: rgb("#fdeaea"),
    inset: (x: 2pt, y: 1.5pt),
    radius: 1pt,
    above: 2.5pt,
    below: 2.5pt,
    [#text(fill: c-eig, weight: "bold")[!] #body],
)
// worked template
#let ex(title, body) = block(
    width: 100%,
    fill: luma(243),
    inset: (x: 2pt, y: 1.8pt),
    radius: 1pt,
    above: 3pt,
    below: 3pt,
    breakable: false,
    [#text(weight: "bold")[#title:] #body],
)
// highlighted sub-block (left rule)
#let focus(c, body) = block(width: 100%, stroke: (left: 1.2pt + c), inset: (left: 2.5pt, y: 1pt), body)

///////////////////
// Formula Sheet //
///////////////////

#align(center)[
    #text(size: 8.5pt, weight: "bold")[ENAE 641 Midterm Formula Sheet] \
    #text(size: 6pt, fill: luma(90))[V. Srivastava · Fall 2026 · HW 01--03 scope]
]

///////////////////////////////
// 1. Linear Algebra By Hand //
///////////////////////////////

#sec([1 #h(0.2em) Linear algebra by hand], c-la)
#let s1 = sub.with(c-la)

#s1[Products.] $A vbu(x) = sum_k x_k vbu(a)_k$ (a combination of the columns $vbu(a)_k$) $= [vbu(alpha)_i^TT vbu(x)]_i$ (dot products with the rows $vbu(alpha)_i^TT$); $(A B)_(i j) = vbu(alpha)_i^TT vbu(b)_j$ (inner sizes match, outer give the size). $(A B)^TT = B^TT A^TT$, $(A B)^(-1) = B^(-1) A^(-1)$, $(A^TT)^(-1) = (A^(-1))^TT$.

#s1[Det & inverse.]
$ mat(a, b; c, d)^(-1) = 1/(a d - b c) mat(d, -b; -c, a), quad A^(-1) = adj(A)/det(A) $
$adj(A)_(i j) = (-1)^(i+j) M_(j i)$ (transposed cofactors). 3×3: expand along any row/col, $det A = sum_j (-1)^(i+j) a_(i j) M_(i j)$ ($M_(i j)$ = minor: delete row $i$, col $j$), or Sarrus for $mat(a, b, c; d, e, f; g, h, i)$: $a e i + b f g + c d h - c e g - a f h - b d i$. $det(A B) = det A det B$, $det A^TT = det A$, $det A^(-1) = 1\/det A$, $det(c A) = c^n det A$; triangular: product of diagonal; $det mat(A, B; 0, D) = det A det D$; swapping rows flips the sign, adding a multiple of a row changes nothing. Blocks: $mat(A, B; 0, D)^(-1) = mat(A^(-1), -A^(-1) B D^(-1); 0, D^(-1))$, $blkdiag{A_i}^(-1) = blkdiag{A_i^(-1)}$.

#s1[Rank, range, nullspace (RREF).] Row-reduce $A$ ($p times m$): rank $rho$ = \# pivots = \# independent columns = \# independent rows ($rho(A^TT) = rho(A)$). $cal(R)(A) subset RR^p$: span of the *pivot columns of the original* $A$. $cal(N)(A) = {vbu(x): A vbu(x) = vbu(0)} subset RR^m$: write pivot variables in terms of free ones; one basis vector per free variable. #key[$rho(A) + nu(A) = m$ (\# columns)]. Square $A$: $det A = 0 <=> rho < n <=> nu > 0 <=> A^(-1)$ DNE $<=>$ columns dependent $<=> lambda = 0$ is an eigenvalue.

#s1[Span, independence, basis.] Stack vectors as columns of $M$: independent $<=> rho(M) =$ \# vectors (square: $det M != 0$); $vbu(y) in spn$ $<=> rho([M, vbu(y)]) = rho(M)$. Any $n$ independent vectors (any nonsingular matrix) form a basis of $RR^n$; coordinates of $vbu(x)$ in basis $V$: solve $V vbu(c) = vbu(x)$, $vbu(c) = V^(-1) vbu(x)$ (unique).

#ex[Template (HW1)][$A = mat(2, 1, 0; 4, 2, 0; 6, 3, 0) ~> mat(1, 1/2, 0; 0, 0, 0; 0, 0, 0)$: $rho = 1$, $nu = 2$, $cal(R) = spn{vec(2, 4, 6)}$, $cal(N) = spn{vec(-1, 2, 0), vec(0, 0, 1)}$.]

#focus(c-la)[
    #s1[Eigenvalues by hand.] $lambda$ with $A vbu(v) = lambda vbu(v)$, $vbu(v) != vbu(0)$ $<=> (lambda I - A)$ singular:
    + Characteristic polynomial $p(lambda) = det(lambda I - A) = 0$ ($n$ roots with multiplicity).
    + 2×2: #key[$lambda^2 - (tr A) lambda + det A = 0$], $lambda = tr/2 plus.minus sqrt((tr\/2)^2 - det)$; complex iff $det > (tr\/2)^2$: $sigma = tr\/2$, $omega = sqrt(det - (tr\/2)^2)$.
    + 3×3: $lambda^3 - (tr A) lambda^2 + (M_(11) + M_(22) + M_(33)) lambda - det A$ ($M_(i i)$: principal 2×2 minors).
    + Roots without a calculator: integer roots divide the constant term (try $plus.minus 1, plus.minus 2, dots$); deflate by synthetic division; then quadratic formula or complete the square $(lambda + a)^2 + b^2 = 0 => lambda = -a plus.minus j b$.
    + Check: $sum lambda_k = tr A$, $product lambda_k = det A$.
    *Shortcuts:* triangular / block-triangular / block-diagonal $=>$ diagonal entries / eigenvalues of the diagonal blocks; companion $mat(0, 1, 0; 0, 0, 1; -a_0, -a_1, -a_2)$ (or its transpose) $=> lambda^3 + a_2 lambda^2 + a_1 lambda + a_0$; $mat(a, b; -b, a) => a plus.minus j b$; $rho(A) = r < n => lambda = 0$ ($>= n - r$ eigvecs); constant row sum $c => lambda = c$, $vbu(v) = [1, dots, 1]^TT$. $A$, $A^TT$, $P^(-1) A P$ share $lambda$; $A^k, A^(-1), A + c I, c A, e^(A t) -> lambda^k, 1\/lambda, lambda + c, c lambda, e^(lambda t)$ (same $vbu(v)$); real $A$: conjugate pairs; symmetric: real $lambda$, orthogonal $vbu(v)$; skew: $j omega$ or 0.

    #s1[Eigenvectors by hand.]
    + For each $lambda$ solve $(A - lambda I) vbu(v) = vbu(0)$ by RREF; one vector per free variable; \# independent $= n - rho(A - lambda I)$ (geometric multiplicity).
    + 2×2 $mat(a, b; c, d)$: $vbu(v) = vec(b, lambda - a)$ or $vec(lambda - d, c)$. 3×3 with $rho(A - lambda I) = 2$: $vbu(v) = vbu(r)_i times vbu(r)_j$ (cross product of two independent rows of $A - lambda I$; fine for complex $lambda$). Companion with 1s above the diagonal: $vbu(v) = [1, lambda, lambda^2]^TT$ (2×2: $[1, lambda]^TT$); for its transpose that vector is the *left* eigenvector.
    + Scale is arbitrary; $overline(lambda) -> overline(vbu(v))$; real pair $Re vbu(v), Im vbu(v)$; check $A vbu(v) = lambda vbu(v)$.
    *Left eigenvectors:* rows $vbu(w)_k^TT$ of $V^(-1)$ (2×2: invert $V$ with the formula above; $vbu(w)_k^TT A = lambda_k vbu(w)_k^TT$, i.e. $A^TT vbu(w)_k = lambda_k vbu(w)_k$), scaled so $vbu(w)_i^TT vbu(v)_j = delta_(i j)$. *Repeated* $lambda$ (alg. mult. $r$): diagonalizable $<=> n - rho(A - lambda I) = r$ for every $lambda$ (always if distinct). Else chain $(A - lambda I) vbu(p)_1 = 0$, $(A - lambda I) vbu(p)_(k+1) = vbu(p)_k$ $=> A P = P J$. \# Jordan blocks $= n - rho(A - lambda I)$; \# blocks of size $>= k$ $= rho((A - lambda I)^(k-1)) - rho((A - lambda I)^k)$.
]

#ex[(a)][$mat(0, 1; -6, -5)$: $lambda^2 + 5 lambda + 6 = (lambda + 2)(lambda + 3)$; $vbu(v) = [1, lambda]^TT$: $vbu(v)_1 = vec(1, -2)$, $vbu(v)_2 = vec(-1, 3)$; $V^(-1) = mat(3, 1; 2, 1)$ $=> vbu(w)_1^TT = [3, 1]$, $vbu(w)_2^TT = [2, 1]$.]
#ex[(b) HW2][$A = mat(0, 0, -13; 1, 0, -17; 0, 1, -5)$: $p = lambda^3 + 5 lambda^2 + 17 lambda + 13$, $p(-1) = 0 => p = (lambda + 1)((lambda + 2)^2 + 9)$, $lambda = -1, -2 plus.minus 3 j$. Rows 2, 3 of $A - lambda I$: $vbu(r)_2 times vbu(r)_3 = [lambda^2 + 5 lambda + 17, lambda + 5, 1]^TT$ $=> vbu(v)_(-2+3j) = [2 + 3j, 3 + 3j, 1]^TT$, $vbu(v)_(-1) = [13, 4, 1]^TT$.]
#ex[(c) repeated][$A = mat(-2, 1; -1, -4)$: $(lambda + 3)^2$; $A + 3 I = mat(1, 1; -1, -1)$, $rho = 1 =>$ one block; $vbu(v)_1 = vec(1, -1)$, $(A + 3I) vbu(p)_2 = vbu(v)_1 => vbu(p)_2 = vec(1, 0)$; $e^(A t) = e^(-3t)[I + t(A + 3 I)]$.]

#s1[Definiteness] (symmetric $M$; only $M_"sym" = (M + M^TT)\/2$ enters $vbu(x)^TT M vbu(x)$). PD $<=>$ all $lambda > 0$ $<=>$ leading principal minors $> 0$ (Sylvester); $mat(a, b; b, c) succ 0 <=> a > 0, a c - b^2 > 0$. PSD $<=>$ all $lambda >= 0$ $<=>$ *all* principal minors $>= 0$.

#s1[Commuting.] $A B = B A$ for: diagonal pairs; $A$ with $I$, $A^k$, $A^(-1)$, $p(A)$, $e^(A t)$; $f(t) M$ with $g(tau) M$. #s1[Calculus.] $dv(, t)(vbu(x)^TT P vbu(x)) = dot(vbu(x))^TT P vbu(x) + vbu(x)^TT dot(P) vbu(x) + vbu(x)^TT P dot(vbu(x))$; $dv(, t) M^(-1) = -M^(-1) dot(M) M^(-1)$; Leibniz $dv(, t) integral_(a(t))^(b(t)) f(t, tau) dd(tau) = f(t, b) dot(b) - f(t, a) dot(a) + integral_a^b pdv(f, t) dd(tau)$.

#s1[Laplace.] $e^(a t) <-> 1/(s - a)$, $t^k e^(a t) <-> k!/(s - a)^(k+1)$, $e^(sigma t) cos omega t <-> (s - sigma)/((s - sigma)^2 + omega^2)$, $e^(sigma t) sin omega t <-> omega/((s - sigma)^2 + omega^2)$, $dot(f) <-> s F - f(0^-)$, $delta <-> 1$; linear only for *constant* coefficients. Cover-up: $c_i = [(s - p_i) F(s)]_(s = p_i)$; double pole: $c_2 = [(s - p)^2 F]_p$, $c_1 = [dv(, s)(s - p)^2 F]_p$; complex pair: write $(alpha (s - sigma) + beta omega)\/((s - sigma)^2 + omega^2)$.

#s1[Integrals.] $integral e^(a t) dd(t) = e^(a t)\/a$; $integral t e^(a t) dd(t) = e^(a t)(a t - 1)\/a^2$; $integral e^(a t) cos b t dd(t) = e^(a t)(a cos b t + b sin b t)\/(a^2 + b^2)$; $integral e^(a t) sin b t dd(t) = e^(a t)(a sin b t - b cos b t)\/(a^2 + b^2)$. #s1[Values.] $sin(pi\/6, pi\/4, pi\/3) = 1\/2, sqrt(2)\/2, sqrt(3)\/2$ (cos reversed); $sqrt(2) approx 1.414$, $sqrt(3) approx 1.732$, $sqrt(5) approx 2.236$; $ln 2 approx 0.693$, $ln 50 approx 3.91$; $e^(-1) approx 0.368$, $e^(-4) approx 0.018$.

////////////////////////////////////////////
// 2. Models, Linearization, and Equilibria //
////////////////////////////////////////////

#sec([2 #h(0.2em) Models, linearization, equilibria], c-mod)
#let s2 = sub.with(c-mod)

#s2[State.] $z_1, dots, z_n$ such that $vbu(z)(t_0)$ plus $vbu(w)(t)$, $t >= t_0$, determine $vbu(eta)(t)$, $t >= t_0$. Implicit $vbu(0) = vbu(f)(dot(vbu(z)), vbu(z), vbu(w), t)$ (explicit $dot(vbu(z)) = vbu(f)(vbu(z), vbu(w), t)$ can be too restrictive), output $vbu(eta) = vbu(h)(vbu(z), vbu(w), t)$. Linear $<=> vbu(f) = E(t) dot(vbu(z)) + F(t) vbu(z) + G(t) vbu(w)$, $vbu(h) = C vbu(z) + D vbu(w)$; standard form ($E$ invertible): #key[$A = -E^(-1) F$, $B = -E^(-1) G$]. Constant matrices $=>$ time-invariant.

#s2[Linearization] about a reference solution $vbu(z)^*(t)$ driven by $vbu(w)^*(t)$: Jacobians
$ E = lr([pdv(vbu(f), dot(vbu(z)))])_*, quad F = lr([pdv(vbu(f), vbu(z))])_*, quad G = lr([pdv(vbu(f), vbu(w))])_*, quad C = lr([pdv(vbu(h), vbu(z))])_* $
$D = [pdv(vbu(h), vbu(w))]_*$, $E_(i j) = pdv(f_i, dot(z)_j)$, etc., evaluated at the reference *after* differentiating. Then
$ dot(vbu(x)) = A(t) vbu(x) + B(t) vbu(u) + vbu(Gamma), quad vbu(y) = C vbu(x) + D vbu(u) + vbu(Gamma)_y $
$vbu(x) = vbu(z) - vbu(z)^*$, $vbu(u) = vbu(w) - vbu(w)^*$, $vbu(y) = vbu(eta) - vbu(eta)^*$; $norm(vbu(Gamma))$ (H.O.T.) small only while $norm(vbu(x)), norm(vbu(u))$ small (Taylor).
#warn[$vbu(x) -> vbu(0)$ means $vbu(z) -> vbu(z)^*$: back to the *reference*, not to 0. $vbu(u)$ is not the physical input: apply $vbu(w) = vbu(w)^* + vbu(u)$; $vbu(y)$ is a deviation too. Reference trajectory $=>$ LTV $A(t)$ (e.g. unicycle); equilibrium $=>$ LTI.]
#s2[Equilibria / trim.] $dot(vbu(z))^* = 0$, $vbu(z)^*$ const: $vbu(0) = vbu(f)(vbu(0), vbu(z)^*, vbu(w)^*)$ is *algebraic*: solve $vbu(z)^*$ given $vbu(w)^*$ or vice versa. One $vbu(w)^*$ can give several equilibria with opposite stability (pendulum $w^* = 0$: $theta^* = 0$ stable, $pi$ unstable).

#ex[Pendulum][$m l^2 dot.double(theta) = -m g l sin theta - b(dot(theta)) + tau_m$, $vbu(z) = [theta, dot(theta)]^TT$, $w = tau_m$:
    $ E = mat(1, 0; 0, m l^2), quad F = mat(0, -1; m g l cos theta^*, b'(dot(theta)^*)), quad G = vec(0, -1) $
    $ => A = mat(0, 1; -alpha_0, -alpha_1), quad B = vec(0, beta_0) $
    $alpha_0 = (g\/l) cos theta^*$, $alpha_1 = b'(dot(theta)^*)\/(m l^2)$, $beta_0 = 1\/(m l^2)$; I/O: $dot.double(y) + alpha_1 dot(y) + alpha_0 y = beta_0 u$. Trim: #key[$m g l sin theta^* = w^*$]. Slide cases ($g\/l = 5$, $b'(0)\/(m l^2) = 2$):
    #table(
        columns: 4,
        align: (center, center, center, left),
        [$theta^*$], [char. poly], [poles], [behaviour],
        [$0$], [$s^2 + 2s + 5$], [$-1 plus.minus 2j$], [decaying oscillation],
        [$pi\/4$], [$s^2 + 2s + 3.54$], [$-1 plus.minus 1.59j$], [same decay, lower $omega$],
        [$80 degree$], [$s^2 + 2s + 0.87$], [$-0.64, -1.36$], [no oscillation],
        [$90 degree$], [$s^2 + 2s$], [$0, -2$], [*not stable*],
        [$180 degree$], [$s^2 + 2s - 5$], [$1.45, -3.45$], [unstable],
    )
]
#image("assets/pendulum-cases.svg")

#s2[Modes & poles.] LTI fundamental solutions $e^(p t) = e^(sigma t)(cos omega t + j sin omega t)$, $p = sigma + j omega$: $abs(e^(p t)) = e^(sigma t) -> 0 <=> Re p < 0$; repeated $p$ $=> t^ell e^(p t)$. 2% settling: #key[$t_s = ln 50 \/ abs(sigma) approx 4\/abs(sigma)$]. $s^2 + 2 zeta omega_n s + omega_n^2$: $p = -zeta omega_n plus.minus j omega_d$, $omega_n = abs(p)$, #key[$zeta = abs(sigma)\/abs(p)$] $= cos(angle "from" -RR)$, $omega_d = omega_n sqrt(1 - zeta^2)$; $zeta < 1$ under-, $= 1$ critically, $> 1$ over-damped. Sketch: envelope $plus.minus A e^(sigma t)$, period $2 pi\/omega_d$, settles to the equilibrium.
#image("assets/splane-modes.svg")
#s2[2×2 type] ($tau = tr A$, $Delta = det A$, $lambda = tau\/2 plus.minus sqrt(tau^2\/4 - Delta)$): AS $<=> tau < 0$ & $Delta > 0$; $Delta < 0$ saddle; $tau = 0 < Delta$ center; $Delta = tau^2\/4$ repeated. $A = mat(0, 1; -k, -b)$ ($dot.double(y) + b dot(y) + k y = 0$): AS $<=> k > 0, b > 0$; $k < 0$ saddle.
#image("assets/trdet.svg")
#ex[PD on a linearization][$dot.double(y) + alpha_1 dot(y) + alpha_0 y = beta_0 u$, $u = -k_p y - k_d dot(y)$ $=> s^2 + (alpha_1 + beta_0 k_d) s + (alpha_0 + beta_0 k_p)$; match $2 zeta omega_n$, $omega_n^2$: $k_d = (2 zeta omega_n - alpha_1)\/beta_0$, $k_p = (omega_n^2 - alpha_0)\/beta_0$. Apply $w = w^* - k_p (theta - theta^*) - k_d dot(theta)$.]
#s2[Hoop] $dot(omega) + b omega + a^2 sin theta = Omega^2 cos theta sin theta$: $A = mat(0, 1; -k, -b)$, $k = a^2 cos theta^* - Omega^2 cos 2 theta^*$; bottom $k = a^2 - Omega^2$, top $k = -a^2 - Omega^2$ (saddle); $Omega^2 > a^2$: bottom becomes a saddle, new $cos theta^* = a^2\/Omega^2$ with $k = Omega^2 sin^2 theta^* > 0$.

//////////////////////////////////////
// 3. Solutions of the State Equation //
//////////////////////////////////////

#sec([3 #h(0.2em) Solutions: fundamental matrix & STM], c-sol)
#let s3 = sub.with(c-sol)

#s3[Problem.] $dot(vbu(x)) = A(t) vbu(x) + B(t) vbu(u)$, $vbu(x)(t_0) = vbu(x)_0$; knowing $n$ solutions from independent ICs gives all solutions.
#s3[Fundamental matrix.] $Psi = [vbu(psi)_1, dots, vbu(psi)_n]$, $dot(vbu(psi))_k = A vbu(psi)_k$, *any* independent ICs: $dot(Psi) = A(t) Psi$, $Psi(t_0)$ nonsingular. $vbu(x)(t) = Psi(t) vbu(alpha)$, $vbu(alpha) = Psi^(-1)(t_0) vbu(x)_0$ conserved when $vbu(u) = vbu(0)$. Liouville: $dv(, t) det Psi = tr A(t) det Psi$ $=> det Psi(t) = det Psi(t_0) e^(integral_(t_0)^t tr A dd(tau)) != 0$ (unless $integral tr A = -oo$).
#s3[STM.] #key[$Phi(t, tau) = Psi(t) Psi^(-1)(tau)$], unique (same for every $Psi$); $vbu(x)(t) = Phi(t, t_0) vbu(x)_0$. Properties: (1) $Phi(t, t) = I$; (2) nonsingular for all $t, tau$; (3) $Phi^(-1)(t, tau) = Phi(tau, t)$; (4) $pdv(, t) Phi(t, tau) = A(t) Phi(t, tau)$, $pdv(, tau) Phi(t, tau) = -Phi(t, tau) A(tau)$ (2nd argument!); (5) $Phi(t_2, t_0) = Phi(t_2, t_1) Phi(t_1, t_0)$ (composition) $=> Phi(t, tau) = Phi(t, 0) Phi(tau, 0)^(-1)$. Also $dv(, t) Psi^(-1) = -Psi^(-1) A$.
#s3[General solution] (zero-input + zero-state):
$ #key[$vbu(x)(t) = underbrace(Phi(t, t_0) vbu(x)_0, vbu(x)_"zi") + underbrace(integral_(t_0)^t Phi(t, tau) B(tau) vbu(u)(tau) dd(tau), vbu(x)_"zs")$] $
$vbu(y) = C(t) Phi(t, t_0) vbu(x)_0 + integral_(t_0)^t G(t, tau) vbu(u)(tau) dd(tau)$, impulse response #key[$G(t, tau) = C(t) Phi(t, tau) B(tau) + D(t) delta(t - tau)$] ($equiv$ using $D(tau) delta(t - tau)$; $p times m$; $integral f(tau) delta(t - tau) dd(tau) = f(t)$). View: $vbu(x) = Psi vbu(alpha)(t)$, $dot(vbu(alpha)) = Psi^(-1) B vbu(u)$.
#s3[Peano--Baker.] $Phi = I + integral_tau^t A(sigma_1) dd(sigma_1) + integral_tau^t A(sigma_1) integral_tau^(sigma_1) A(sigma_2) dd(sigma_2) dd(sigma_1) + dots$
#key[*If* $A(t) M = M A(t)$, $M(t, tau) = integral_tau^t A(sigma) dd(sigma)$, *then* $Phi(t, tau) = e^(M(t, tau)) = sum_k M^k\/k!$]
#warn[Only when $A(t)$ commutes with its integral; in general $Phi != exp integral A$.]
#s3[Commuting sum] $A(t) = sum_j f_j(t) A_j$, $A_i A_j = A_j A_i$: $Phi = product_j e^(A_j integral_tau^t f_j (sigma) dd(sigma))$. #s3[Block-triangular] $mat(A_1, A_(12); 0, A_2) => Phi = mat(Phi_1, Phi_(12); 0, Phi_2)$, $Phi_(12)(t, tau) = integral_tau^t Phi_1(t, sigma) A_(12)(sigma) Phi_2(sigma, tau) dd(sigma)$. #s3[Scalar:] $dot(x) = a(t) x + b(t) u => x = e^(integral_(t_0)^t a) x_0 + integral_(t_0)^t e^(integral_sigma^t a) b(sigma) u(sigma) dd(sigma)$. #s3[Coordinates] $vbu(x) = P(t) vbu(z)$: $Phi_x = P(t) Phi_z (t, tau) P^(-1)(tau)$.
#s3[LTI.] $Phi = e^(A(t - tau))$ (elapsed time only $=> t_0 = 0$): #key[$vbu(x) = e^(A t) vbu(x)_0 + integral_0^t e^(A(t - tau)) B vbu(u)(tau) dd(tau)$], $G(t) = C e^(A t) B + D delta(t)$.
#s3[Building $Phi$ for LTV by hand:] (i) triangular/decoupled rows: solve the scalar row, substitute upward; (ii) check commutativity $=>$ exponential of the integral; (iii) split $A = sum f_j A_j$; (iv) find $P(t)$ making $A_p$ constant (e.g. a rotation); (v) periodic $A(t)$: Floquet; or any $n$ independent solutions $=> Phi = Psi(t) Psi^(-1)(tau)$.
#ex[Triangular (HW2)][$dot(x)_1 = t x_2$, $dot(x)_2 = -x_2$: $x_2 = e^(-(t - tau)) x_2(tau)$, $x_1 = x_1(tau) + integral_tau^t sigma e^(-(sigma - tau)) dd(sigma) thin x_2(tau)$ $=> Phi = mat(1, (tau + 1) - (t + 1) e^(-(t - tau)); 0, e^(-(t - tau)))$.]
#ex[Commuting (HW2)][$mat(cos 3t, sin 3t; -sin 3t, cos 3t) = cos 3t thin I + sin 3t thin S$, $S = mat(0, 1; -1, 0)$, $S^2 = -I$: $Phi = e^((sin 3t - sin 3tau)\/3) mat(cos theta, sin theta; -sin theta, cos theta)$, $theta = (cos 3tau - cos 3t)\/3$.]

///////////////////////////////////////////
// 4. Matrix Exponential, Resolvent, TFM //
///////////////////////////////////////////

#sec([4 #h(0.2em) Matrix exponential, resolvent, TFM], c-exp)
#let s4 = sub.with(c-exp)

#s4[Definition.] $e^(A t) = sum_(k=0)^oo t^k/k! A^k = Phi(t, 0)$. $e^(A 0) = I$; nonsingular; $e^(A(t_1 + t_2)) = e^(A t_1) e^(A t_2)$; $(e^(A t))^(-1) = e^(-A t)$; $(e^(A t))^TT = e^(A^TT t)$; $dv(, t) e^(A t) = A e^(A t) = e^(A t) A$; $e^(P^(-1) A P t) = P^(-1) e^(A t) P$; eigenvalues $e^(lambda_k t)$.
#warn[$e^M != [e^(m_(i j))]$ (not entrywise). $e^((A + B) t) = e^(A t) e^(B t)$ *only if* $A B = B A$.]
#table(
    columns: (auto, 1fr),
    align: (center + horizon, left + horizon),
    [$A$], [$e^(A t)$],
    [$diag{a_k}$], [$diag{e^(a_k t)}$],
    [$mat(a, b; -b, a)$], [$e^(a t) mat(cos b t, sin b t; -sin b t, cos b t)$],
    [$mat(a, 1; 0, a)$], [$e^(a t) mat(1, t; 0, 1)$],
    [$J_3 (a)$], [$e^(a t) mat(1, t, t^2\/2; 0, 1, t; 0, 0, 1)$; $J_n$: $t^k\/k!$ on superdiag. $k$],
    [$blkdiag{A_i}$], [$blkdiag{e^(A_i t)}$],
    [$theta S$, $S^2 = -I$], [$cos theta thin I + sin theta thin S$],
    [$N$, $N^k = 0$], [$I + N t + dots + N^(k-1) t^(k-1)\/(k-1)!$ (finite)],
)
#s4[2×2 closed forms.] Distinct real $lambda_(1,2)$: $e^(A t) = (lambda_1 e^(lambda_2 t) - lambda_2 e^(lambda_1 t))/(lambda_1 - lambda_2) I + (e^(lambda_1 t) - e^(lambda_2 t))/(lambda_1 - lambda_2) A$; repeated $lambda$: $e^(lambda t)[I + t(A - lambda I)]$; complex $sigma plus.minus j omega$: #key[$e^(sigma t)[cos omega t thin I + (sin omega t)/omega (A - sigma I)]$].
#s4[Methods.]
+ *Resolvent:* #key[$e^(A t) = cal(L)^(-1){(s I - A)^(-1)}$], $(s I - A)^(-1) = adj(s I - A)\/det(s I - A) = Q(s)\/r(s)$: every entry has the *same* denominator $r(s)$ (degree $n$; numerators degree $<= n - 1$) $=>$ partial fractions $=>$ every entry is a combination of the same $e^(lambda_k t)$: poles = eigenvalues; repeated $=> t^ell e^(lambda t)$.
+ *Cayley--Hamilton:* $p_A (A) = 0 => A^(k >= n)$ is a combination of $I, dots, A^(n-1)$ $=> e^(A t) = sum_(k=0)^(n-1) alpha_k (t) A^k$, $sum_k alpha_k lambda_i^k = e^(lambda_i t)$ (Vandermonde, *$lambda$ distinct*; repeated: add $dv(, lambda)$: $sum_k k alpha_k lambda^(k-1) = t e^(lambda t)$). Mostly a theory tool.
+ *Diagonalize:* $A^k = V Lambda^k V^(-1) =>$ #key[$e^(A t) = V e^(Lambda t) V^(-1) = sum_k e^(lambda_k t) vbu(v)_k vbu(w)_k^TT$] (dyads $vbu(v)_k vbu(w)_k^TT$, $n times n$).
+ *Jordan / block modal:* $e^(A t) = P e^(J t) P^(-1)$, block by block.
#ex[Two ways][$A = mat(0, 1; -6, -5)$: $(s I - A)^(-1) = 1/((s + 2)(s + 3)) mat(s + 5, 1; -6, s)$, e.g. $(s + 5)/((s + 2)(s + 3)) = 3/(s + 2) - 2/(s + 3)$. Dyads (eigs from §1a): $vbu(v)_1 vbu(w)_1^TT = mat(3, 1; -6, -2)$, $vbu(v)_2 vbu(w)_2^TT = mat(-2, -1; 6, 3)$ $=> e^(A t) = e^(-2t) mat(3, 1; -6, -2) + e^(-3t) mat(-2, -1; 6, 3)$.]
#s4[TFM.] $(s I - A) vbu(X) = vbu(x)_0 + B vbu(U)$, $vbu(Y) = C(s I - A)^(-1) vbu(x)_0 + G(s) vbu(U)$:
$ #key[$G(s) = C(s I - A)^(-1) B + D$] quad (p times m) $
$G(t) = cal(L)^(-1){G(s)} = C e^(A t) B + D delta(t)$; $[G(s)]_(i j) = Y_i\/U_j$ (zero IC, other inputs 0); modal partial fractions #key[$[G(s)]_(i j) = sum_(k=1)^n ((vbu(c)_i^TT vbu(v)_k)(vbu(w)_k^TT vbu(b)_j))/(s - lambda_k) + d_(i j)$]; invariant under $vbu(x) = P vbu(z)$.
#warn[Laplace is LTI only: $cal(L){A(t) vbu(x)(t)} != A(t) vbu(X)(s)$ (nor $cal(L){A} vbu(X)$).]

////////////////////////////////////////////
// 5. Eigenstructure, Modal & Coordinates //
////////////////////////////////////////////

#sec([5 #h(0.2em) Modal forms, coordinates, Jordan, Floquet], c-eig)
#let s5 = sub.with(c-eig)

#s5[Eigenstructure.] $A V = V Lambda$; distinct $lambda$ $=>$ independent $vbu(v)_k$ $=> V^(-1)$ exists, $A = V Lambda V^(-1)$. $V^(-1) = mat(vbu(w)_1^TT; dots.v; vbu(w)_n^TT)$: left eigenvectors ($A^TT$ has the same $lambda$); #key[$vbu(w)_i^TT vbu(v)_j = delta_(i j)$] (bi-orthonormal).
#s5[Modal response.] #key[$vbu(x)(t) = sum_k (vbu(w)_k^TT vbu(x)_0) e^(lambda_k t) vbu(v)_k$]: modes $e^(lambda_k t) vbu(v)_k$, weights $vbu(w)_k^TT vbu(x)_0$ = how much the IC excites mode $k$. $vbu(x)_0 = beta vbu(v)_j => vbu(x) = beta e^(lambda_j t) vbu(v)_j$: only mode $j$, straight along the $vbu(v)_j$ line.
#s5[Forced, decoupled:] $vbu(x) = sum_k z_k vbu(v)_k = V vbu(z)$ with #key[$dot(z)_k = lambda_k z_k + sum_j (vbu(w)_k^TT vbu(b)_j) u_j$], $y_i = sum_k (vbu(c)_i^TT vbu(v)_k) z_k + vbu(d)_i^TT vbu(u)$. $vbu(w)_k^TT vbu(b)_j$: strength of input $j$ on mode $k$ ($vbu(w)_k perp vbu(b)_j$: input $j$ cannot excite mode $k$); $vbu(c)_i^TT vbu(v)_k$: visibility of mode $k$ in output $i$ ($perp$: invisible). Modal form: $A_m = V^(-1) A V = Lambda$, $B_m = V^(-1) B$ (row $k$ = $vbu(w)_k^TT B$), $C_m = C V$, $D_m = D$; same I/O behaviour.
#image("assets/slide-modal-coords.png", width: 62%)
#s5[Real block modal form.] $lambda = sigma plus.minus j omega$ with $vbu(v)$ for $sigma + j omega$: $A [Re vbu(v), Im vbu(v)] = [Re vbu(v), Im vbu(v)] mat(sigma, omega; -omega, sigma)$; $P = [Re vbu(v)_1, Im vbu(v)_1, Re vbu(v)_3, Im vbu(v)_3, vbu(v)_5, dots]$ $=> A_p = blkdiag(mat(sigma_1, omega_1; -omega_1, sigma_1), dots, lambda_5, dots)$, $e^(A_p t)$ = rotation blocks.
*One oscillatory mode:* $vbu(x)_0 = a Re vbu(v) + b Im vbu(v)$ $=> z_1 = e^(sigma t)(a cos omega t + b sin omega t)$, $z_2 = e^(sigma t)(b cos omega t - a sin omega t)$, $vbu(x) = z_1 Re vbu(v) + z_2 Im vbu(v)$: only that frequency, and $vbu(x)(t) in spn{Re vbu(v), Im vbu(v)}$ for all $t$ (a plane). Real mode: $vbu(x)_0 prop vbu(v)$. General IC: $vbu(z)_0 = P^(-1) vbu(x)_0$. Polar $v_i = r_i e^(j phi_i)$: $x_i$ oscillates with relative amplitude $r_i$ and phase $phi_i$. Which input excites a mode: that mode's rows of $B_p = P^(-1) B$.
#ex[HW2 (eigs from §1b)][$P = [Re vbu(v), Im vbu(v), vbu(v)_3] = mat(2, 3, 13; 3, 3, 4; 1, 0, 1)$, $A_p = P^(-1) A P = mat(-2, 3, 0; -3, -2, 0; 0, 0, -1)$; $vbu(x)_0 = [2, 3, 1]^TT$ excites only the oscillation, $[13, 4, 1]^TT$ only $e^(-t)$.]
#s5[Coordinate change] $vbu(x) = P vbu(z)$: $A_p = P^(-1) A P$, $B_p = P^(-1) B$, $C_p = C P$, $D_p = D$; time-varying #key[$A_p (t) = P^(-1)(A P - dot(P))$] (most general). Invariant: $Psi_x = P Psi_z$, $Phi_x (t, tau) = P(t) Phi_z (t, tau) P^(-1)(tau)$, $G(t, tau)$, $G(s)$ (insert $I = P P^(-1)$), eigenvalues (constant $P$).
#s5[Jordan form.] $det(s I - A) = (s - lambda_1)^r dots$: $P^(-1) A P = blkdiag{J_1(lambda_1), dots, J_ell (lambda_1), lambda_(r+1), dots}$, $J_i = mat(lambda, 1, ; , dots.down, 1; , , lambda)_(n_i times n_i)$, $sum n_i = r$ (algebraic), #key[$ell = nu(lambda_1 I - A) = n - rho(lambda_1 I - A)$] (geometric = \# eigvecs = \# blocks). $ell = r$: diagonal; $ell = 1$: one $r times r$ block; $1 < ell < r$: sizes not fixed ($r = 4$: $ell = 3 => (2,1,1)$; $ell = 2 => (3,1)$ or $(2,2)$), set by chains $(lambda I - A)^k vbu(w) = 0$. Largest block $n_i$ $=>$ terms up to $t^(n_i - 1) e^(lambda t)$. Very sensitive: $mat(1, 1; epsilon, 1)$ is one block at $epsilon = 0$, diagonal for $epsilon != 0$.
#ex[HW2 Jordan][$(lambda^2 + 2 lambda + 10)(lambda + 3)^3$, $rho(A + 3I) = 4 => ell = 1$, one 3×3 block: $P^(-1) e^(A t) P = blkdiag(e^(-t) mat(cos 3t, sin 3t; -sin 3t, cos 3t), e^(-3t) mat(1, t, t^2\/2; 0, 1, t; 0, 0, 1))$; $vbu(x)_0 = vbu(p)_4 => P^(-1) vbu(x) = [0, 0, t e^(-3t), e^(-3t), 0]^TT$.]
#s5[Lyapunov transformation.] Any $A(t)$ $->$ any constant $R$ via $P(t) = Psi(t) e^(-R t)$ ($A_p = R$), so $Phi = P(t) e^(R(t - tau)) P^(-1)(tau)$. But $dot(x) = x -> dot(z) = -z$ with $P = e^(2t)$ ("unstable made stable": $P^(-1) -> 0$). Lyapunov transformation: also $abs(det P(t)) >= eta > 0$ and $abs(det P^(-1)(t)) >= eta$ $forall t$ $=>$ preserves convergence / divergence.
#s5[Floquet] ($A(t + T) = A(t)$): monodromy #key[$Q = Phi(t_0 + T, t_0)$], $R = 1/T log Q$ ($e^(R T) = Q$), $P(t) = Phi(t, t_0) e^(-R t)$ is $T$-periodic (a Lyapunov transformation) since $Phi(t + T, tau + T) = Phi(t, tau)$:
$ #key[$Phi(t, tau) = P(t) e^(R(t - tau)) P^(-1)(tau)$] quad "(periodic part" times "LTI part)" $
AS $<=> Re lambda(R) < 0 <=> abs(mu) < 1$ for all multipliers $mu = lambda(Q)$; any $abs(mu) > 1$: unstable. $log Q = V diag(ln mu_k) V^(-1)$; $mu < 0 => ln mu = ln abs(mu) + j pi$ (complex $R$; or use period $2T$).
#warn[Frozen-time eigenvalues of $A(t)$ do *not* decide LTV stability: HW2's periodic $A(t)$ has $lambda = -1\/4 plus.minus j sqrt(7)\/4$ for all $t$, yet $vbu(x) = e^(t\/2)[cos t, -sin t]^TT$.]

///////////////////////////
// 6. Stability, Lyapunov //
///////////////////////////

#sec([6 #h(0.2em) Stability & Lyapunov], c-stab)
#let s6 = sub.with(c-stab)

#s6[Definitions] ($vbu(x) = Phi(t, t_0) vbu(x)_0$).
- *Stable i.s.L.:* $exists gamma(t_0) < oo$: $norm(vbu(x)(t)) <= gamma norm(vbu(x)_0)$ $forall t >= t_0$ $<=> norm(Phi(t, t_0))_i$ bounded; *uniform* if $gamma$ is independent of $t_0$.
- *AS:* i.s.L. *and* $norm(vbu(x)) -> 0$ ($<=> norm(Phi) -> 0$): $forall epsilon$ $exists T(epsilon, t_0)$: $norm(vbu(x)(t)) <= epsilon norm(vbu(x)_0)$ $forall t >= t_0 + T$. i.s.L. is needed: a finite escape time could blow up and then decay. *UAS:* $T$ independent of $t_0$ (still depends on $epsilon$).
- *ES:* #key[$norm(vbu(x)(t)) <= gamma e^(-lambda(t - t_0)) norm(vbu(x)_0)$], $lambda > 0$ ($norm(vbu(x))$ may wiggle under the envelope); *UES:* $gamma, lambda$ independent of $t_0$.
UES $=>$ ES $=>$ AS $=>$ i.s.L.; AS $arrow.r.not$ ES. #key[Linear: UAS $<=>$ UES]: $norm(Phi(t_0 + T, t_0)) <= epsilon$ $forall t_0$ $=> norm(vbu(x)) <= gamma epsilon^k norm(vbu(x)_0)$ on $[t_0 + k T, t_0 + (k + 1) T]$ $=> lambda = -ln epsilon \/ T$.
#image("assets/slide-uas-ues.png")
#s6[Counterexamples.] $dot(x) = (4 t sin t - 2t) x$: $abs(x(t_0 + pi))\/abs(x_0) = e^((2 t_0 + pi)(4 - pi)) -> oo$ as $t_0 -> oo$: $gamma(t_0)$ unbounded, stability not uniform (slide writes $sin 2t$; its ratio needs $sin t$). $Phi = (1 + t_0)\/(1 + t)$: AS but not UAS ($epsilon = 1\/2$: $T = 1 + t_0$, so $T(1\/2, 0) = 1$, $T(1\/2, 10) = 11$), not ES (algebraic decay).
#s6[LTI.] AS $<=>$ ES $<=>$ all $Re lambda < 0$ (Hurwitz); i.s.L. $<=>$ all $Re lambda <= 0$ and $j omega$-axis eigenvalues have 1×1 Jordan blocks (else $t^k$ growth). LTV: use $Phi$, Lyapunov, or Floquet (never frozen eigenvalues).
#s6[Norms.] $norm(vbu(x))_1 = sum abs(x_k)$, $norm(vbu(x))_2 = sqrt(sum x_k^2)$, $norm(vbu(x))_oo = max abs(x_k)$, $norm(vbu(x))_p = (sum abs(x_k)^p)^(1\/p)$; $norm(vbu(x)) >= 0$, $= 0 <=> vbu(x) = 0$, $norm(alpha vbu(x)) = abs(alpha) norm(vbu(x))$, triangle. Equivalence: $underline(c) norm(vbu(x))_q <= norm(vbu(x))_p <= overline(c) norm(vbu(x))_q$ (stability is norm-independent). Induced $norm(A)_i = max norm(A vbu(x))\/norm(vbu(x))$, #key[$norm(A vbu(x)) <= norm(A)_i norm(vbu(x))$]: $norm(A)_(i,2) = sqrt(lambda_max (A^TT A))$, $norm(A)_(i,1)$ = max abs. column sum, $norm(A)_(i,oo)$ = max abs. row sum; $norm(A + B) <= norm(A) + norm(B)$, $norm(A B) <= norm(A) norm(B)$.
#s6[Quadratic forms.] $vbu(x)^TT M vbu(x) = vbu(x)^TT M_"sym" vbu(x)$ since $vbu(x)^TT M_"skew" vbu(x) = 0$; Rayleigh--Ritz (symmetric $M$): #key[$lambda_min (M) norm(vbu(x))_2^2 <= vbu(x)^TT M vbu(x) <= lambda_max (M) norm(vbu(x))_2^2$]. $M(t)$ PD (PSD): $vbu(x)^TT M vbu(x) > 0$ ($>= 0$) $forall t$; uniformly PD: $min_t lambda_min (M(t)) >= eta > 0$.
#s6[First test] ($P = I$): $dv(, t) norm(vbu(x))_2^2 = vbu(x)^TT (A^TT + A) vbu(x) <= alpha(t) norm(vbu(x))_2^2$, $alpha = lambda_max (A^TT + A)$; #key[$max_t lambda_max (A^TT (t) + A(t)) <= 0 =>$ i.s.L.] ($gamma = 1$).
#warn[Sufficient only: $lambda_max > 0$ does *not* imply instability (too restrictive: $norm(vbu(x))$ may rise and still decay). Skew $A(t)$ ($A^TT = -A$) $=> norm(vbu(x))_2$ constant: uniformly i.s.L., not AS (HW3).]
#s6[Lyapunov functions.] Nested closed level sets $V(vbu(x)) = c$; $dot(V) = pdv(V, vbu(x)) dot dot(vbu(x)) <= 0$ (outward normal $dot$ velocity) $=>$ motion points inward $=> S(c) = {vbu(x): V <= c}$ is *invariant* $=>$ i.s.L. Any nested closed contours work; linear systems: ellipses, QLF $V = vbu(x)^TT P vbu(x)$, $P succ 0$ ($P = I$ is $norm(vbu(x))^2$).
#image("assets/lyap-levelsets.svg")
#s6[Lyapunov equation.] $P(t)$ bounded, uniformly PD: $rho = max_t lambda_max (P)$, $eta = min_t lambda_min (P)$:
$ #key[$-Q(t) = dot(P) + A^TT P + P A$] quad => quad dot(V) = -vbu(x)^TT Q vbu(x) $
*Theorem:* $Q succ.eq 0 =>$ i.s.L. ($eta norm(vbu(x))^2 <= V(t) <= V(t_0) <= rho norm(vbu(x)_0)^2$, $gamma = sqrt(rho\/eta)$). $Q$ uniformly PD, $nu = min_t lambda_min (Q) > 0$ $=>$ ES: $dot(V) <= -nu norm(vbu(x))^2 <= -(nu\/rho) V$ and, by comparison,
$ #key[$norm(vbu(x)(t))_2 <= sqrt(rho/eta) e^(-nu/(2 rho) (t - t_0)) norm(vbu(x)_0)_2$] $
Failed test $=>$ *inconclusive* (try another $P$). Comparison lemma: $dot(v) <= f(v, t)$, $v(t_0) <= w(t_0)$, $dot(w) = f(w, t)$ $=> v(t) <= w(t)$.
#s6[Converse.] $A(t)$ UES $=>$ for *any* bounded unif. PD $Q$ there is a bounded unif. PD $P(t) = integral_t^oo Phi^TT (tau, t) Q(tau) Phi(tau, t) dd(tau)$ ($norm(P) <= gamma_Q c^2 \/ (2 lambda)$); one pair $=>$ infinitely many.
#s6[LTI:] #key[$A^TT P + P A = -Q$]: $A$ Hurwitz $=>$ $forall Q succ 0$ $exists P succ 0$.
#warn[Does *not* say every $P succ 0$ gives $Q succ 0$: $P = I$ (the $lambda(A + A^TT)$ test) can fail for a stable non-normal $A$ (HW3 F16); infinitely many other $P$ work.]
#ex[2×2 ALE][$A = mat(0, 1; -k, -b)$, $Q = I$, $P = mat(p_1, p_2; p_2, p_3)$: (1,1) $-2 k p_2 = -1$; (1,2) $p_1 - b p_2 - k p_3 = 0$; (2,2) $2(p_2 - b p_3) = -1$ $=> p_2 = 1/(2k)$, $p_3 = (1 + k)/(2 b k)$, $p_1 = (1 + k)/(2b) + b/(2k)$; $P succ 0 <=> k, b > 0$; then $rho = lambda_max (P)$, $eta = lambda_min (P)$, $nu = 1$.]
#s6[Sign-changing $nu(t)$] (HW3): $dot(V) <= -nu(t) norm(vbu(x))^2 => dot(V) <= -beta V$, $beta = nu^+\/rho - nu^-\/eta$, $nu^(plus.minus) = max{plus.minus nu, 0}$ (lecture: $V <= V_0 e^(-1/rho integral alpha)$, valid where $alpha >= 0$) $=> norm(vbu(x)) <= sqrt(rho\/eta) e^(-1/2 integral_(t_0)^t beta) norm(vbu(x)_0)$. $integral beta$ bounded below (e.g. $integral nu^- < oo$) $=>$ i.s.L.; $integral beta -> oo =>$ AS; $integral_t^(t + T) beta >= mu > 0$ $forall t => $ UES. $nu = 1\/t$: AS only; $-c e^(-t)$: i.s.L.; $1\/t^2$: i.s.L. only ($integral nu < oo$); $2 - 3 sin t$: UES if $rho\/eta lt.tilde 12.3$. $dot(V) = alpha V$, $alpha > 0$, gives $V -> oo$ only if $integral alpha = oo$.
#s6[Local stability.] $dot(vbu(x)) = A(t) vbu(x) + vbu(Gamma)(vbu(x), t)$, $vbu(Gamma)$ 2nd order: $forall epsilon$ $exists delta(epsilon)$: $norm(vbu(Gamma)) <= epsilon norm(vbu(x))$ on $B_delta = {norm(vbu(x)) <= delta}$. With $(P, Q)$ of the linear part, $gamma_p = max_t norm(P)_i$:
$ dot(V) <= -nu norm(vbu(x))^2 + 2 gamma_p norm(vbu(x)) norm(vbu(Gamma)) <= -(nu - 2 gamma_p epsilon) norm(vbu(x))^2 = -nu_1 norm(vbu(x))^2 $
in $B_delta$ if $epsilon < nu\/(2 gamma_p)$. Pick $c$ with $S_c = {V <= c} subset B_delta$ (e.g. $c = eta delta^2$; then $S_c supset B_(delta sqrt(eta\/rho))$): $S_c$ invariant, $vbu(x) -> 0$ exponentially from any $vbu(x)_0 in S_c$. *Local*, not global; $S_c$ is a conservative estimate of the region of attraction; needs a UES linearization. Linearization with some $Re lambda > 0$: unstable; with $Re lambda = 0$: inconclusive. Physically $vbu(z) -> vbu(z)^*$.
#image("assets/slide-roa.png", width: 50%)
#ex[Which tool][
    - *LTI stable?* eigenvalues ($Re lambda < 0$), or the ALE with any $Q succ 0$.
    - *LTV stable?* not frozen eigenvalues: bound $norm(Phi)$, a Lyapunov $P(t)$, or Floquet multipliers (periodic).
    - *$e^(A t)$:* 2×2 closed forms; dyads $V e^(Lambda t) V^(-1)$; resolvent + partial fractions.
    - *One mode:* $vbu(x)_0$ in its eigen-subspace; input $j$ reaches mode $k$ iff $vbu(w)_k^TT vbu(b)_j != 0$.
    - *Nonlinear equilibrium:* trim $->$ Jacobians $->$ $lambda(A)$: local only, $vbu(z) -> vbu(z)^*$.
]
