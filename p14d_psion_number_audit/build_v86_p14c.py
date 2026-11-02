p='/home/claude/maia-karlosmarden/psion-proofs/p14d_psion_number_audit/paper_v86/PH_Foundations_v86.tex'
src=open(p).read()
assert 'P14-C' not in src, 'ja aplicado'
def rep(old,new,count=1):
    global src
    n=src.count(old); assert n==count,(n,old[:80]); src=src.replace(old,new)
rep("%      tables gain P14-B; abstract and summaries updated.\n",
"%      tables gain P14-B; abstract and summaries updated.\n%   3. App. L.4.7: hadron-size geometry recorded as phenomenological; P14-C opened\n%      (baryon binding scale and internal structure consistent with J.Sigma mixing).\n")
row_b=r"P14-B & Conserved baryon (rung) label; sign and size of the bifurcation bias (§13.4, §13.7) & open & $\eta$B, nucleon stability, Theorem KS \\ \hline"
row_c=r"P14-C & Baryon binding scale ($(\hslash c)/a\approx163$ MeV) and internal structure consistent with $J\!\cdot\!\Sigma$ mixing; hexagonal radii are phenomenological (App.~L.4.7) & open & hadron radii, $m_n-m_p$ \\ \hline"
rep(row_b, row_b+"\n"+row_c, count=2)
l47=r"""

\noindent \textbf{L.4.7 Hadron size and geometry (v86, P14-C, [T4]).} An exploratory, pre-registered series (companion repository, WP-GEO to WP-GEO12 and WP-LJS) placed the six labels on geometric arrangements and computed charge-weighted mean-square radii with the charges $Q_k$ and equal constituent masses. The planar Z6 ring in natural label order (the hexagon, and equivalently any arrangement that keeps the cyclic order on one circle) is the only arrangement tested that reproduces, with one scale fixed by the proton, $R_n=\langle r^2\rangle_n/\langle r^2\rangle_p=-0.154$ (measured $-0.164$) and $R_{\Sigma^-}=0.846$ ($0.86\pm0.21$), and it gives $\langle r^2\rangle_{\Sigma^0}=0$ exactly; the scale is $a=1.212$ fm. Of the 720 orderings of the labels on the ring, only the natural cyclic order passes. The same arrangement fails for the triples containing $k=3$ ($\Lambda$, $\Sigma^+$, against chiral-perturbation-theory and lattice predictions), fails for magnetic moments, and cannot make leptons pointlike: leptonic sizes must come from the contact scale $1/\Lambda_{\rm meta}$, consistent with the exact $J\!\cdot\!\Sigma$-inertness of $(0,1,2)$. Three results prevent promotion. (a) $J\!\cdot\!\Sigma$ conserves $\sum k$ and therefore mixes $n$, $\Sigma^0$, $\Lambda$ (and $p$ with the two $\Sigma^+$ candidates); at the coupling of this paper the mixing destroys $R_n$ with adiabatic labels, and with ground-state labels gives $R_n\approx-0.08$ and $m_n<m_p$. (b) The scale $(\hslash c)/a=163$ MeV is not derived: in four dimensions the $J\!\cdot\!\Sigma$ contact term would need tuning to one part in $10^{24}$, and none of four pre-registered transmutation forms built from the constants of this paper lands within a factor 2 (nearest: rung $\ell=8$, factor 3.7; $e^{-S/g^2}$, factor 14). (c) The $k=3$ position that fits $\Lambda$ and $\Sigma^+$ (a $30^\circ$ shift on a dodecagon) is a discrete fit with no dynamical origin. This is recorded as open problem \textbf{P14-C}: a derived baryon binding scale and an internal structure consistent with the $J\!\cdot\!\Sigma$ dynamics. Nothing in the main text depends on it."""
i=src.index(r"\noindent \textbf{L.4.6 Dual sector and nucleon stability (v86).}")
j=src.index("\n\n",i)
src=src[:j]+l47+src[j:]
rep("open-problem tables gain P14-B;", "open-problem tables gain P14-B and P14-C; App.~L.4.7 added (hadron-size geometry, phenomenological);")
open(p,'w').write(src); print("ok")
