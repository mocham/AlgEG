# Depends_on

1. [completely reducible residual representations are quasi-semisimple](F09-012.md)

# Info


* Unit class: Theorem
* Status: Verified
* Management: Explicit
* Sage verification: None.
* Title: Rank-one connected-center groups have crystalline lifts without central descent
* Scope: main-spine
* Scope: rank-one
* Scope: crystalline-lifts
* Keyword: connected-center
* Keyword: pgl2
* Keyword: gl2
* Keyword: central-descent-bypass
* Source: [Lin22 quasi-semisimplicity and crystalline lifting, lines 866-924 and 1281-1460](Provenance/Lin22/source/G-irr-journal.tex)
* Source: [EG23 strong existence of crystalline lifts, lines 26618-26704](Provenance/EG23/moduli.tex)

# Statement


Let $F/\mathbb Q_p$ be finite with $p>2$, and let $\mathbf G$ be a split connected reductive group scheme with connected geometric center and semisimple rank at most one. Every continuous representation $\bar\rho:G_F\to\mathbf G(\overline{\mathbb F}_p)$ admits a crystalline lift $\rho:G_F\to\mathbf G(\mathcal O_{\overline{\mathbb Q}_p})$.

# Proof or Examples


## Expanded self-contained proof

Let $X=X^*(\mathbf T)$ be the character lattice of a split maximal torus. If the semisimple rank is zero, $\mathbf G$ is a split torus. Suppose it is one, with roots $\{\pm\alpha\}$. Connectedness of the center says that $X/\mathbb Z\alpha$ is torsion-free, so $\alpha$ is primitive. Choose a basis $(\alpha,e_2,\ldots,e_r)$ of $X$. The coroot has coordinates $\alpha^\vee=(2,b_2,\ldots,b_r)$. Replacing the $e_i$ by integral combinations and by $e_i+n_i\alpha$ reduces $(b_2,\ldots,b_r)$ either to zero or to $(1,0,\ldots,0)$, according as all $b_i$ are even or not. The two resulting root data are those of $\operatorname{PGL}_2\times\mathbb G_m^{r-1}$ and $\operatorname{GL}_2\times\mathbb G_m^{r-2}$, respectively. Thus it suffices to treat $\mathbb G_m$, $\operatorname{GL}_2$, and $\operatorname{PGL}_2$.

The strong crystalline-lift theorem of EG23, Theorem 6.4.4, gives a crystalline lift of every residual $\operatorname{GL}_d$ representation, hence handles $\mathbb G_m$ and $\operatorname{GL}_2$. For $\operatorname{PGL}_2$, first suppose that $\bar\rho$ is contained in a Borel. The normalization of an upper-triangular projective matrix by its lower-right entry identifies that Borel with the subgroup $\{\left(\begin{smallmatrix}t&u\\0&1\end{smallmatrix}\right)\}\subset\operatorname{GL}_2$. Apply the EG23 theorem to the resulting residual $\operatorname{GL}_2$ representation and projectivize its crystalline lift. If $\bar\rho$ is not contained in a Borel, it is $\operatorname{PGL}_2$-irreducible and therefore completely reducible. F09-012 makes it quasi-semisimple, and Lin22, Theorem 1281-1301, gives a crystalline $\operatorname{PGL}_2$ lift. Finally, take the Cartesian product of the lifts under the displayed direct-product decomposition of $\mathbf G$. No lift through a finite central cover and no vanishing assertion in $H^2(G_F,Z)$ is used.

## Faithful GrobV14 proof

GrobV14 does not print this fact as a separately labeled unit. Its managed role is the result titled "Rank-one connected-center groups have crystalline lifts without central descent". It is a direct proof dependency of `COR18-002`, `F05-027`.

## Post-proof remarks

No internal issues. The connected-center rank-one root-datum reduction is explicit and has only the PGL2-times-torus and GL2-times-torus cases. EG23 lines 26618-26704 prove the required GL_d pointwise theorem. Lin22 lines 866-924 and 1281-1460 prove quasi-semisimplicity and crystalline lifting for the irreducible PGL2 case. The Borel-contained PGL2 case uses the explicit algebraic upper-triangular GL2 section before applying EG23. F09-012 is Verified. This route bypasses, rather than assumes away, the central-isogeny obstruction.
