# Depends_on

None.

# Info


* Unit class: Theorem
* Status: Verified
* Management: Explicit
* Sage verification: None.
* Title: Regular crystalline lifts for GL2, PGL2, and their Borels
* Slogan: Every residual GL2 or PGL2 point has a regular crystalline lift, and a Borel-valued point has one in the same Borel.
* Scope: main-spine
* Scope: rank-one
* Scope: crystalline-lifts
* Keyword: GL2
* Keyword: PGL2
* Keyword: Borel
* Keyword: regular-Hodge-type
* Keyword: trivial-inertial-type
* Source: [EG23 extension theorem, ordinary definition, and strong regular crystalline lifting theorem, lines 25731-25776,26546-26562,26618-26682](Provenance/EG23/moduli.tex)
* Source: [Lin22 regularity definition and regular crystalline lifting theorem, lines 1531-1553,1689-1767](Provenance/Lin22/source/G-irr-journal.tex)

# Statement


Let $K/\mathbb Q_p$ be finite and let $G$ be $\operatorname{GL}_2$ or $\operatorname{PGL}_2$. Every continuous representation $\bar\rho:G_K\to G(\overline{\mathbb F}_p)$ has a crystalline lift $\rho:G_K\to G(\mathcal O_{\overline{\mathbb Q}_p})$ of regular labeled Hodge type and trivial inertial type. If $\bar\rho$ takes values in a Borel subgroup $B$ of the special fiber, then $\rho$ can be chosen to take values in a Borel of the split integral model reducing to $B$.

# Proof or Examples


## Expanded self-contained proof

For $G=\operatorname{GL}_2$, EG23 Theorem 6.4.4 gives a crystalline lift of regular labeled Hodge--Tate weights. If the residual representation is contained in a Borel, use its invariant line to write it as the specified extension of one character by another. The rank-two instance of the induction in the proof of Theorem 6.4.4 applies EG23 Theorem 6.4.3 to that chosen extension and produces an integral extension of crystalline characters. Thus the chosen line, not merely some line after generic-fiber conjugacy, lifts. Equivalently, the result is ordinary in the exact sense of EG23 Definition 6.4.2 and is upper triangular in a basis lifting the prescribed residual basis. Its labeled weights are regular by Theorem 6.4.4.

Let $q:\operatorname{GL}_2\to\operatorname{PGL}_2$ be the quotient and let $B_{\operatorname{PGL}_2}$ be the standard Borel. The subgroup $H=\{\left(\begin{smallmatrix}t&u\\0&1\end{smallmatrix}\right):t\in\mathbb G_m, u\in\mathbb G_a\}$ of $\operatorname{GL}_2$ maps isomorphically to $B_{\operatorname{PGL}_2}$. The inverse is the integral group-scheme morphism that normalizes an upper-triangular projective matrix by its lower-right entry. Hence a residual $B_{\operatorname{PGL}_2}$-valued representation has the genuine homomorphic lift $s\circ\bar\rho$ to $H$. Apply the Borel-preserving $\operatorname{GL}_2$ case and projectivize. Crystallinity and inertial type are functorial under $q$. If the two labeled $\operatorname{GL}_2$ weights are distinct, the unique root of $\operatorname{PGL}_2$ evaluates on the projected Hodge cocharacter as their nonzero difference, so regularity is also preserved. The projectivized lift remains in $B_{\operatorname{PGL}_2}$.

If a residual $\operatorname{PGL}_2$ representation is not contained in a Borel, it is $\operatorname{PGL}_2$-irreducible because the Borels are the only proper parabolics. It is therefore $G$-completely reducible, and Lin22 Theorem 1690--1701 gives a Hodge--Tate regular crystalline lift. Lin22 identifies this regularity with regularity of the labeled Hodge cocharacters. Finally, any nonstandard residual Borel can be conjugated to the standard one; smoothness lifts the conjugating point to the coefficient ring, so conjugating back gives a Borel reducing to the original $B$. Every crystalline lift has trivial Weil--Deligne inertia, which gives the asserted inertial type in all cases.

## Faithful GrobV14 proof

GrobV14 does not print this fact as a separately labeled unit. Its managed role is the result titled "Regular crystalline lifts for GL2, PGL2, and their Borels". Its direct managed consumers are `D04-005`, `D08-001`; none is one of the 78 indexed in-paper units.

## Post-proof remarks

Verified exact EG23 regular GL2 and chosen-extension Borel lifting, the integral PGL2-Borel section and projectivization, and Lin22 regular lifting for the non-Borel irreducible case.
